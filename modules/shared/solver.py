import scikits.odes.sundials.ida as ida
from numpy import (zeros, concatenate, all, sum, power, isscalar, linspace,
                   asarray, cumsum, empty)

simulation_err_str = ('%s simulation: Calculation did not reach '
                      'the expected time. An error occured.'
                      'Reached time: % 10.6f\nExpected time: % 10.6f')


def simulate_direct(initialize_z0, model, residual_fn,
                    update_initial_condition=None,
                    root_fn = None, nr_rootfns=None, algvars_idx=None):

    # Check supplied arguments
    if (not root_fn is None) and (nr_rootfns is None):
        raise Exception("Error: Function 'root_fn' was set in "
                        "simulate_inverse, but not the number of rooting "
                        "functions 'nr_rootfns'.")

    # Initialization
    verbosity = model.verbosity

    measurements_nr = model.iterations + 1

    t   = empty([measurements_nr, ], dtype=float)
    z   = empty([measurements_nr, model.z_size], dtype=float)

    t0 = 0.0

    t[0]    = t0

    if update_initial_condition:
        z0 = empty([model.z_size, ], float)
    else:
        z0 = z[0, :]
    zp0 = zeros(z0.shape, float)

    i = 1
    model.set_iteration(i)   # set model for the first iteration
    initialize_z0(z0, model) # and compute the initial state

    if update_initial_condition: # as initial state can be modified at each
        z[0, :] = z0             # restart, z0 was preallocated separately

    solver = ida.IDA(residual_fn,
                     compute_initcond='yp0',
                     first_step_size=model.first_step_size,
                     atol=model.atol, rtol=model.rtol,
                     max_step_size = model.max_step_size,
                     max_steps = model.max_steps,
                     algebraic_vars_idx=algvars_idx,
                     rootfn=root_fn, nr_rootfns=nr_rootfns,
                     #linsolver='band', uband=1, lband=1,
                     user_data=model)

    solver_initialized = False

    # Run computation
    out_s = '{: >5d}. {: 12.1f}  {: 12.1f}  {: 12.1f}'

    iterations = model.iterations

    if verbosity > 1:
        capt_s = '{:>6} {:>12}  {:>12}  {:>12}'
        print(capt_s.format('Run', 'Start time', 'Duration', 'End time'))
    while True:
        if verbosity == 2:
            total_duration = (model.duration + model.fh_duration
                              + model.deceleration_duration)
            print(out_s.format(i, t0, total_duration, t0 + total_duration))

        for duration_type in ['duration', 'fh_duration',
                              'deceleration_duration']:

            duration = getattr(model, duration_type)

            if duration == 0.0: continue

            t_end = t0 + duration

            solver.set_options(tstop=t_end)

            if verbosity > 2:
                print(out_s.format(i, t0, duration, t_end))

            if duration_type == 'duration':
                model.set_omega2g_fn('centrifugation')
            elif duration_type == 'fh_duration':
                model.set_omega2g_fn('falling_head')
                # backup values
                r0           = model.r0
                duration     = model.duration

                model.r0    = model.r0_fall
                model.duration = fh_duration
            else:
                model.set_omega2g_fn('deceleration')

            if not solver_initialized:
                solver.init_step(t0, z0, zp0)
                solver_initialized = True

            (flag, t_out) = solver.step(t_end, z[i, :])

            if t_out < t_end:
                if verbosity > 1:
                    print('Error occured during computation. Solver failed '
                          'at time\nt_err=', t_out,
                          '\nExpected value of t:', t)
                if verbosity > 2:
                    print('Values at this time:\nz_err=', z_retn)

                return (False, t[:i], z[:i, :])


            t0 = t_out

            if duration_type == 'fh_duration': # restore backuped values
                model.r0    = r0
                model.duration = duration

        t[i] = t_out

        if i == iterations: break

        # update values for the next iteration
        i = i+1
        model.set_iteration(i)

        if update_initial_condition:
            z0[:] = z[i-1, :]
            update_initial_condition(i, z0, model)
            if not all(z0 == z[i-1, :]):
                solver_initialized = False
        else:
            z0 = z[i-1, :]

        if model.always_restart_solver:
            solver_initialized = False

    return (True, t, z)

def simulate_inverse(times, direct_fn, model, init_parameters,
                     wl_in_meas=None, wl_out_meas=None,
                     gc_meas=None, rm_meas=None,
                     wl_in_weights=None, wl_out_weights=None,
                     gc_weights=None, rm_weights=None,
                     optimfn='leastsq'):

    from modules.shared.functions import determine_scaling_factor
    from modules.shared.show import disp_inv_results
    from numpy import log, exp, alen

    available_solvers = ['leastsq', 'fmin', 'fmin_power', 'fmin_cg',
                         'fmin_bfgs', 'raster']
    if not optimfn in available_solvers:
        print("Unknown inverse method solver for 'optimfn': ", optimfn)
        print("Available solvers are: ", available_solvers)
        exit(1)


    transform = {'ks': lambda ks: log(ks),
                 'n':  lambda n: log(n - 1.0),
                 'gamma': lambda gamma: log(-gamma)}
    untransform = {'ks': lambda ks_transf: exp(ks_transf),
                   'n': lambda n_transf: 1+exp(n_transf),
                   'gamma': lambda gamma_transf: -exp(gamma_transf)}

    optimized_parameters = []

    def update_model(optim_params, model):
        for (name, value) in zip(optimized_parameters, optim_params):
            setattr(model, name, untransform[name](value))

        if 'n' in optimized_parameters:
            model.m = 1 - 1/model.n

        if model.verbosity > 0:
            print()
            for name in optimized_parameters:
                if name == 'ks':
                    print('{:5}: {: .8g}'.format('Ks', getattr(model, name)))
                else:
                    print('{:5}: {: .8g}'.format(name, getattr(model, name)))

    lbounds = {}
    ubounds = {}

    def penalize(model, when='out_of_bounds'):
        penalization = 0.0

        if when == 'out_of_bounds':
            for param in optimized_parameters:
                value = getattr(model, param)
                if lbounds[param] > value:
                    a = exp(value - lbounds[param])
                    penalization = (penalization + 10 * (a + 1/a))
                elif ubounds[param] < value:
                    a = exp(value - ubounds[param])
                    penalization = (penalization + 10 * (a + 1/a))
        else:
            for param in optimized_parameters:
                value = getattr(model, param)
                a = exp(value - lbounds[param])
                b = exp(value - ubounds[param])

                penalization = (penalization + 10 * (a + 1/a) + 10 * (b + 1/b))

        return penalization

    calc_wl_in  = bool(wl_in_meas)
    calc_wl_out = bool(wl_out_meas)
    calc_gc     = bool(gc_meas)
    calc_rm     = bool(rm_meas)

    if not ((wl_in_weights is None) and (wl_out_weights is None)
            and (gc_weights is None) and (rm_weights is None)):
        add_weights = True
    else:
        add_weights = False

    no_measurements = empty([0,], dtype=float)

    if calc_wl_in:
        wl_in_M = asarray(wl_in_meas, dtype=float)
        wl_in_scale_coef = determine_scaling_factor(wl_in_meas)
        wl_in_M[:] = wl_in_M * wl_in_scale_coef
    else:
        wl_in_M = no_measurements

    if calc_wl_out:
        wl_out_M = cumsum(asarray(wl_out_meas, dtype=float))
        wl_out_scale_coef = determine_scaling_factor(wl_out_M)
        wl_out_M[:] = wl_out_M * wl_out_scale_coef
    else:
        wl_out_M = no_measurements

    if calc_gc:
        gc_M = np.asarray(gc_meas, dtype=float)
        gc_scale_coef = determine_scaling_factor(gc_meas)
        gc_M[:] = gc_M * gc_scale_coef
    else:
        gc_M = no_measurements

    if calc_rm:
        rm_M = np.asarray(rm_meas, dtype=float)
        rm_scale_coef = determine_scaling_factor(rm_meas)
        rm_M[:] = rm_M * rm_scale_coef
    else:
        rm_M = no_measurements

    measurements = concatenate((wl_in_M, wl_out_M, gc_M, rm_M))

    if add_weights:
        if wl_in_weights is None:
            wl_in_wghts = no_measurements
        else:
            wl_in_wghts = asarray(wl_in_weights)

        if wl_out_weights is None:
            wl_out_wghts = no_measurements
        else:
            wl_out_wghts = asarray(wl_out_weights)

        if gc_weights is None:
            gc_wghts = no_measurements
        else:
            gc_wghts = asarray(gc_weights)

        if rm_weights is None:
            rm_wghts = no_measurements
        else:
            rm_wghts = asarray(rm_weights)

        weights = concatenate((wl_in_wghts, wl_out_wghts, gc_wghts, rm_wghts))

    iteration = 0

    def optimfn_wrapper(optimargs):
        nonlocal iteration

        print(15 * '*', ' Iteration: {:4d}'.format(iteration), ' ', 15 * '*')
        iteration += 1

        update_model(optimargs, model)

        penalization = penalize(model, when='out_of_bounds')

        if penalization > 0.0:
            if model.verbosity > 1:
                print('Optimized arguments are out of bounds... Penalizing by ',
                      penalization)

            if optimfn == 'leastsq':
                return penalization + measurements
            else:
                return penalization * alen(measurements)

        (flag, t, wl_in, wl_out, gc, rm) = direct_fn(model)

        if flag:
            # direct computation went O.K.
            if calc_wl_in:
                wl_in_C = wl_in * wl_in_scale_coef
            else:
                wl_in_C = no_measurements

            if calc_wl_out:
                wl_out_C = wl_out * wl_out_scale_coef
            else:
                wl_out_C = no_measurements

            if calc_gc:
                gc_C = gc * gc_scale_coef
            else:
                gc_C = no_measurements

            if calc_rm:
                rm_C = rm * rm_scale_coef
            else:
                rm_C = no_measurements

            if model.verbosity > 0:
                disp_inv_results(model, t, inv_params=None,
                                 wl_in_inv=wl_in, wl_out_inv=wl_out,
                                 gc1_inv=gc, rm1_inv=rm,
                                 display_graphs=False,
                                 disp_abserror=(model.verbosity > 1))

            computation = concatenate((wl_in_C, wl_out_C, gc_C, rm_C))

            error = computation - measurements

            if add_weights:
                error[:] = error * weights

            if optimfn == 'leastsq':
                return error
            else:
                return sum(power(error, 2))

        else:
            # something is wrong, so penalize
            penalization = min(penalize(model, when='always'), 1e10)
            if model.verbosity > 1:
                print('Direct problem did not converge for given optimization '
                      'parameters... Penalizing by ', penalization)

            if optimfn == 'leastsq':
                return penalization + measurements
            else:
                return penalization * alen(measurements)

            return penalization

    init_values = []

    for (param, value) in init_parameters.items():
        if not value is None:
            optimized_parameters.append(param)
            (init_value, (lbound, ubound)) = value

            init_values.append(transform[param](init_value))
            lbounds[param] = lbound
            ubounds[param] = ubound

    import scipy.optimize

    optimize = getattr(scipy.optimize, optimfn)

    if optimfn == 'leastsq':
        inv_params, cov = optimize(optimfn_wrapper, init_values,
                                   epsfcn=model.epsfcn, factor=model.factor,
                                   xtol=model.xtol, ftol=model.ftol)
    else:
        inv_params, cov = optimize(optimfn_wrapper, init_values,
                                   xtol=model.xtol, ftol=model.ftol)

    # we now assume, that in the last run were used the optimal parameters
    # and therefore are still set in the model
    optim_params = {name: getattr(model, name) for name in optimized_parameters}
    model_verbosity = model.verbosity # backup verbosity
    model.verbosity = 0
    (flag, t, wl_in, wl_out, gc, rm) = direct_fn(model)
    model.verbosity = model_verbosity # restore verbosity
    disp_inv_results(model, t, inv_params=optim_params, cov=cov,
                     wl_in_inv=wl_in, wl_out_inv=wl_out, gc1_inv=gc, rm1_inv=rm,
                     disp_abserror=True)

    return inv_params

def simulate_inverse_old(direct_fn, xdata, ydata, init_params,
                     optimfn='leastsq'):

    def lsq_wrapper_fn(xdata, *optim_args):
        direct_results = direct_fn(xdata, optim_args)
        return concatenate(direct_results)

    def leastsq_wrapper_fn(optim_args, *xdata):
        direct_results = direct_fn(xdata[0], optim_args)
        return (concatenate(direct_results) - ydata)

    def fmin_wrapper_fn(optim_args, *xdata):
        direct_results = direct_fn(xdata[0], optim_args)
        tmp = concatenate(direct_results)
        return sum(power(tmp - ydata, 2))

    if optimfn == 'lsq':
        from scipy.optimize import curve_fit
        return curve_fit(lsq_wrapper_fn, xdata, ydata, p0 = init_params,
                         epsfcn=1e-5,
                         factor=0.1)
    if optimfn == 'leastsq':
        from scipy.optimize import leastsq
        return leastsq(leastsq_wrapper_fn, init_params, args=(xdata,),
                       epsfcn=1e-5,
                       factor=0.1)
    elif optimfn == 'fmin':
        from scipy.optimize import fmin
        return fmin(fmin_wrapper_fn, init_params, args=(xdata,))
    elif optimfn == 'fmin_powell':
        from scipy.optimize import fmin_powell
        return fmin_powell(fmin_wrapper_fn, init_params, args=(xdata,))
    elif optimfn == 'fmin_cg':
        from scipy.optimize import fmin_cg
        return fmin_cg(fmin_wrapper_fn, init_params, args=(xdata,))
    elif optimfn == 'fmin_bfgs':
        from scipy.optimize import fmin_bfgs
        return fmin_bfgs(fmin_wrapper_fn, init_params, args=(xdata,))
    elif optimfn == 'raster':
        lbounds = xdata.inv_lbounds
        ubounds = xdata.inv_ubounds
        raster_size = xdata.raster_grid_size
        nrlevels = len(lbounds)
        if isscalar(raster_size):
            raster_size = [raster_size] * nrlevels

        optimargs = [0]* nrlevels
        grid = []
        output = []

        from copy import copy
        for (lb, ub, rsize) in zip(lbounds, ubounds, raster_size):
            grid.append(linspace(float(lb), float(ub), float(rsize)))

        def loopsingle(level):
            for (ind, val) in enumerate(grid[level]):
                optimargs[level] = val

                if level == nrlevels - 1:
                    output.append((fmin_wrapper_fn(optimargs, xdata),
                                  copy(optimargs)))
                else:
                    loopsingle(level+1)

        loopsingle(0)

        print(output)
        from  pickle import dump
        dump(output, '/home/archetyp/Desktop/dump.bin')
