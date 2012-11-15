from __future__ import print_function
import scikits.odes.sundials.ida as ida
from numpy import (zeros, concatenate, all, sum, power, isscalar, linspace,
                   asarray, cumsum, empty, ones)
from modules.shared.functions import has_data

simulation_err_str = ('%s simulation: Calculation did not reach '
                      'the expected time. An error occured.'
                      'Reached time: % 10.6f\nExpected time: % 10.6f')


def simulate_direct(initialize_z0, model, residual_fn,
                    update_initial_condition=None, initialize_zp0=None,
                    root_fn = None, nr_rootfns=None, algvars_idx=None,
                    on_phase_change = None, continue_on_root_found=None):

    # Check supplied arguments
    if (not root_fn is None) and (nr_rootfns is None):
        raise Exception("Error: Function 'root_fn' was set in "
                        "simulate_inverse, but not the number of rooting "
                        "functions 'nr_rootfns'.")

    # Initialization
    verbosity = model.verbosity

    measurements_times = model.measurements_times

    measurements_nr = len(measurements_times)

    t   = empty([measurements_nr, ], dtype=float)
    z   = empty([measurements_nr, model.z_size], dtype=float)

    t0 = 0.0

    model.set_value('omega_start', 0.0)
    model.set_value('t0', t0)
    model.set_value('phase', 'U') # initialize to dummy value

    z0 = empty([model.z_size, ], float)
    zp0 = zeros(z0.shape, float)

    initialize_z0(z0, model) # compute the initial state
    if not initialize_zp0 is None: initialize_zp0(zp0, z0, model)

    if measurements_times[0] == t0:
        t[0]    = t0
        z[0, :] = z0
        i = 1
    else:
        i = 0
    t_meas = measurements_times[i]

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
    previous_phase  = None

    # Run computation
    out_s = '{: >5d}. {: 12.1f}  {: 12.1f}  {: 12.1f}  {:>6}'

    if verbosity > 1:
        capt_s = '\n{:>6} {:>12}  {:>12}  {:>12}  {:>6}'
        print(capt_s.format('Run', 'Start time', 'Duration', 'End time',
                            'Phases'))
    while True:
        if verbosity == 2:
            total_duration = (model.duration + model.fh_duration
                              + model.deceleration_duration)
            phases = []
            if model.duration > 0.0: phases.append('A')
            if model.deceleration_duration > 0.0: phases.append('D')
            if model.fh_duration > 0.0: phases.append('G')
            print(out_s.format(i, t0, total_duration, t0 + total_duration,
                               '-'.join(phases)))

        for phase in ('a', 'd', 'g'):
            if phase == 'a':
                duration = getattr(model, 'duration')
            elif phase == 'd':
                duration = getattr(model, 'deceleration_duration')
            else:
                duration = getattr(model, 'fh_duration')

            if duration == 0.0: continue

            model.t0 = t0
            t_end = t0 + duration

            solver.set_options(tstop=t_end)

            model.phase = phase
            if phase == 'a':
                pass
            elif phase == 'd':
                pass
            else: # phase == 'g'
                # backup values
                backup = model.get_parameters(['re', 'duration'])

                (model.r0, model.duration) = (model.r0_fall, duration)

            model.set_omega2g_fn(phase)

            if verbosity > 2:
                print(out_s.format(i, t0, duration, t_end, phase.upper()))

            if not solver_initialized:
                (flag, t_init) = solver.init_step(t0, z0, zp0)
                if not flag:
                    return (False, t, z)
                solver_initialized = True

            if not on_phase_change is None:
                if previous_phase is None: previous_phase = phase
                if previous_phase != phase: on_phase_change(model, phase)

            while True:
                (flag, t_out) = solver.step(t_meas, z[i, :])

                if flag < 0:     # error occured
                    if verbosity > 1:
                        print('Error occured during computation. Solver failed '
                              'at time\nt_err=', t_out,
                              '\nExpected value of t:', t)
                    if verbosity > 2:
                        print('Values at this time:\nz_err=', z[i, :])

                    return (False, t[:i], z[:i, :])
                elif flag == 2: # root found
                    if ((not continue_on_root_found is None)
                         and (continue_on_root_found(model, t_out, z[i, :]))):
                        if verbosity > 1:
                             print('Root found: continuing computation')
                    else:
                        if verbosity > 1:
                             print('Root found: aborted further computation.')
                        return (False, t[:i], z[:i, :])
                else: # success or t_stop reached
                    if (flag == 0) or (t_out == t_meas):
                        t[i] = t_out
                        i += 1
                        if i < measurements_nr:
                            t_meas = measurements_times[i]
                        else:
                            # previous measurement was the last measurement we
                            # have taken so we can abort further computation
                            break

                    # Assuming t_out is t_end+numerical_error
                    if (flag == 1) or (t_end == t_out):
                        break

            t0 = t_out

            if phase == 'a':
                model.omega_start = model.omega
            elif phase == 'd':
                model.omega_start = 0.0
            else:
                # restore backuped values for 'fh_duration'
                model.set_parameters(backup)

        if not model.next_iteration(): break

        # update values for the next iteration
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

def simulate_inverse(direct_fn, model, init_parameters,
                     measurements, measurements_weights={},
                     optimfn='leastsq'):

    from modules.shared.show import display_status, mk_status_item
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

    add_weights = False
    for weight in measurements_weights.values():
        if not weight is None:
            add_weights = True
            break

    measurements_names = measurements.get_names()
    data_M = measuremetns.get_values()
    measurements_scales = \
       measurements.get_scales(scaling_coefs=model.measurements_scale_coefs)
    (measurements_names, data_M, ) = ([], [], [])
    weights = measurements.get_weights()

    data_scale_coef = concatenate(measurements_scales)
    measurements_sc = concatenate(data_M) * data_scale_coef

    if weights:
        add_weights = True
        weights = concatenate(weights)
    else:
        add_weights = False

    global ITERATION # Python 2.7 hack (no support for nonlocal variables)
    ITERATION = 0

    def optimfn_wrapper(optimargs):
        global ITERATION # Python 2.7 hack (no support for nonlocal variables)

        print(15 * '*', ' Iteration: {:4d}'.format(ITERATION), ' ', 15 * '*')
        ITERATION += 1

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
    iters  = ITERATION

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
    out = ''
    for (name, value) in results:
        if not value is None:
            out += ' |{:>8}'.format(name)
    out += ' |'
    for (name, value) in results:
        if not value is None:
            out += ' |{:8d}'.format(value)
    print(out, '|')

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
