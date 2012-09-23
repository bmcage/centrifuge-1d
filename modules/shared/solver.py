import scikits.odes.sundials.ida as ida
from numpy import (zeros, concatenate, all, sum, power, isscalar, linspace,
                   asarray, cumsum, empty, ones)

simulation_err_str = ('%s simulation: Calculation did not reach '
                      'the expected time. An error occured.'
                      'Reached time: % 10.6f\nExpected time: % 10.6f')


def simulate_direct(initialize_z0, model, residual_fn,
                    update_initial_condition=None, initialize_zp0=None,
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

    model.set_value('omega_start', 0.0)
    model.set_value('t0', t0)

    if update_initial_condition:
        z0 = empty([model.z_size, ], float)
    else:
        z0 = z[0, :]

    initialize_z0(z0, model) # compute the initial state

    zp0 = zeros(z0.shape, float)
    if not initialize_zp0 is None: initialize_zp0(zp0, z0, model)

    i = 1
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
    model.init_iteration() # Re-initialize iterable variables

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

        for duration_type in ['duration', 'deceleration_duration',
                              'fh_duration']:

            duration = getattr(model, duration_type)

            if duration == 0.0: continue

            model.t0 = t0
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
                (flag, t_init) = solver.init_step(t0, z0, zp0)
                if not flag:
                    return (False, t, z)
                solver_initialized = True

            (flag, t_out) = solver.step(t_end, z[i, :])

            if t_out < t_end:
                if verbosity > 1:
                    print('Error occured during computation. Solver failed '
                          'at time\nt_err=', t_out,
                          '\nExpected value of t:', t)
                if verbosity > 2:
                    print('Values at this time:\nz_err=', z[i, :])

                return (False, t[:i], z[:i, :])


            t0 = t_out

            if duration_type == 'duration':
                model.omega_start = model.omega
            elif duration_type == 'deceleration':
                model.omega_start = 0.0
            else:
                # restore backuped values for 'fh_duration'
                model.r0    = r0
                model.duration = duration

        t[i] = t_out

        if not model.next_iteration(): break

        # update values for the next iteration
        i = i+1

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


def set_optimized_variables(optim_params, model, untransform=None):
    """ optim_params: dict of pairs param_name: param_value """
    if untransform:
        for (name, value) in optim_params.items():
            setattr(model, name, untransform[name](value))
    else:
        for (name, value) in optim_params.items():
            setattr(model, name, value)

    if 'n' in optim_params:
        model.m = 1 - 1/model.n

    if model.verbosity > 0:
        print()
        for name in optim_params.keys():
            if name == 'ks':
                print('{:5}: {: .8g}'.format('Ks', getattr(model, name)))
            else:
                print('{:5}: {: .8g}'.format(name, getattr(model, name)))

def simulate_inverse(times, direct_fn, model, init_parameters,
                     wl_in_meas=None, wl_out_meas=None,
                     gc_meas=None, rm_meas=None,
                     wl_in_weights=None, wl_out_weights=None,
                     gc_weights=None, rm_weights=None,
                     optimfn='leastsq'):

    from modules.shared.functions import determine_scaling_factor
    from modules.shared.show import disp_status as display_status, \
     mk_status_item
    from numpy import log, exp, alen

    available_solvers = ['leastsq', 'fmin', 'fmin_powell', 'fmin_cg',
                         'fmin_bfgs', 'raster']
    if not optimfn in available_solvers:
        print("Unknown inverse method solver for 'optimfn': ", optimfn)
        print("Available solvers are: ", available_solvers)
        exit(1)

    max_value = 1e150

    transform = {'ks': lambda ks: max(log(ks), -max_value),
                 'n':  lambda n: max(log(n - 1.0), -max_value),
                 'gamma': lambda gamma: max(log(-gamma), -max_value)}
    untransform = {'ks': lambda ks_transf: min(exp(ks_transf), max_value),
                   'n': lambda n_transf: 1+min(exp(n_transf), max_value),
                   'gamma': lambda gamma_transf: -min(exp(gamma_transf), max_value)}

    optimized_parameters = []

    lbounds = {}
    ubounds = {}

    def penalize(model, when='out_of_bounds'):
        max_penalization = 1e50

        penalization = 0.0

        if when == 'out_of_bounds':
            for param in optimized_parameters:
                value = getattr(model, param)
                if lbounds[param] > value:
                    a = min(exp(value - lbounds[param]), max_penalization)
                    penalization = (penalization
                                    + min(10 * (a + 1/a), max_penalization))
                elif ubounds[param] < value:
                    a = min(exp(value - ubounds[param]), max_penalization)
                    penalization = (penalization
                                    + min(10 * (a + 1/a), max_penalization))
        else:
            for param in optimized_parameters:
                value = getattr(model, param)
                a = min(exp(value - lbounds[param]), max_penalization)
                b = min(exp(value - ubounds[param]), max_penalization)

                penalization = \
                  (penalization + min(10 * (a + 1/a) + 10 * (b + 1/b),
                                      max_penalization))

        return penalization

    if not ((wl_in_weights is None) and (wl_out_weights is None)
            and (gc_weights is None) and (rm_weights is None)):
        add_weights = True
    else:
        add_weights = False

    (measurements_names, data_M, measurements_scales) = ([], [], [])
    measurements_weights = []
    for (data, name, weights) in zip((wl_in_meas, wl_out_meas, gc_meas, rm_meas),
                                     ('MI', 'MO', 'GC', 'RM'),
                                     (wl_in_weights, wl_out_weights, gc_weights,
                                      rm_weights)):
        if not bool(data): continue

        measurement = asarray(data, dtype=float)
        if name == 'MO':
            measurement = cumsum(measurement)
        data_scale_coef = determine_scaling_factor(measurement)
        data_scale = data_scale_coef * ones(measurement.shape, dtype=float)

        measurements_names.append(name)
        data_M.append(measurement)
        measurements_scales.append(data_scale)

        if add_weights:
            if not weights is None:
                measurements_weights.append(asarray(weights, dtype=float))
            else:
                measurements_weights.append(ones(measurement.share, dtype=float))

    data_scale_coef = concatenate(measurements_scales)
    measurements_sc = concatenate(data_M) * data_scale_coef

    if add_weights:
        weights = concatenate(measurements_weights)

    iteration = 0

    def optimfn_wrapper(optimargs):
        nonlocal iteration

        print(15 * '*', ' Iteration: {:4d}'.format(iteration), ' ', 15 * '*')
        iteration += 1

        set_optimized_variables(dict(zip(optimized_parameters, optimargs)),
                                model, untransform=untransform)

        penalization = penalize(model, when='out_of_bounds')

        if penalization > 0.0:
            if model.verbosity > 1:
                print('Optimized arguments are out of bounds... Penalizing by ',
                      penalization)

            if optimfn == 'leastsq':
                return penalization + measurements_sc
            else:
                return penalization * alen(measurements_sc)

        direct_data = direct_fn(model, measurements_names)
        (flag, t, data_C) = (direct_data[0], direct_data[1], direct_data[2:])

        if flag:
            # direct computation went O.K.
            if model.verbosity > 0:
                status_items = []
                for (name, value_C, value_M) in zip(measurements_names,
                                                    data_C, data_M):
                    status_items.append(mk_status_item(name, value_C, value_M))

                display_status(data_plots=status_items)

            computation_sc = data_scale_coef * concatenate(data_C)
            error = computation_sc - measurements_sc

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
                return penalization + measurements_sc
            else:
                return penalization * alen(measurements_sc)

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

    # Initialize output variables that are not present for every solver
    msg = None
    cov = None
    gopt = None
    gcalls = None
    fcalls = None
    iters  = iteration

    # Run optimization
    if optimfn == 'leastsq':
        (opt_params, cov, infodic, msg, ier) = \
          optimize(optimfn_wrapper, init_values,
                   epsfcn=model.epsfcn, factor=model.factor,
                   xtol=model.xtol, ftol=model.ftol,
                   full_output=True)
        fcalls = infodic['nfev']
    elif optimfn == 'fmin':
        (opt_params, fopt, iters, fcalls, warnflag) = \
          optimize(optimfn_wrapper, init_values,
                   xtol=model.xtol, ftol=model.ftol, maxfun=model.max_fev,
                   maxiter=model.max_inv_iter, disp=model.disp_inv_conv,
                   full_output=True, retall=False)
    elif optimfn == 'fmin_powell':
        (opt_params, fopt, direc, iters, fcalls, warnflag) = \
          optimize(optimfn_wrapper, init_values,
                   xtol=model.xtol, ftol=model.ftol, maxfun=model.max_fev,
                   maxiter=model.max_inv_iter, disp=model.disp_inv_conv,
                   full_output=True, retall=False)
    elif optimfn == 'fmin_cg':
        (opt_params, fopt, fcalls, gcalls, warnflag) = \
          optimize(optimfn_wrapper, init_values,
                   maxiter=model.max_inv_iter, gtol=model.gtol,
                   disp=model.disp_inv_conv,
                   full_output=True, retall=False)
    elif optimfn == 'fmin_bfgs':
        (opt_params, fopt, gopt, Bopt, fcalls, gcalls, warnflag) = \
          optimize(optimfn_wrapper, init_values,
                   maxiter=model.max_inv_iter, gtol=model.gtol,
                   disp=model.disp_inv_conv,
                   full_output=True, retall=False)


    # Display inverse solver statistic
    print('\nInverse problem statistics:\n')
    if not msg is None:
        print('\n', msg)
    if not gopt is None:
        print('\nGradient at optimum:\n', gopt, '\n')

    results = [('iters', iters), ('fcalls', fcalls), ('gcalls', gcalls)]
    for (name, value) in results:
        if not value is None:
            print(' |{:>8}'.format(name), end='')
    print(' |')
    for (name, value) in results:
        if not value is None:
            print(' |{:8d}'.format(value), end='')
    print(' |')

    # Run experiment once more with optimal values to display results at optimum
    set_optimized_variables(dict(zip(optimized_parameters, opt_params)),
                            model, untransform)
    optim_params = {name: getattr(model, name) for name in optimized_parameters}

    return (optim_params, cov)

def compute_raster(model):
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
