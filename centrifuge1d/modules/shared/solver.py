from __future__ import print_function, division
import scikits.odes.sundials.ida as ida
from scikits.odes.sundials.ida import StatusEnumIDA
import numpy as np
import sys
from .show import show_inbetween_result
from matplotlib.cbook import flatten

simulation_err_str = ('%s simulation: Calculation did not reach '
                      'the expected time. An error occured.'
                      'Reached time: % 10.6f\nExpected time: % 10.6f')


def simulate_direct(initialize_z0, model, measurements, residual_fn,
                    update_initial_condition=None, initialize_zp0=None,
                    root_fn = None, nr_rootfns=None, algvars_idx=None,
                    exclude_algvar_from_error=False,
                    on_phase_change = None, continue_on_root_found=None,
                    on_measurement=None, truncate_results_on_stop=True,
                    show_progress=True):
    """
    Input:
      initialize_z0 - function with signature f(z0, model), where:
                        z0 - vector to be filled with initial values
                        model - argument supplied to simulate_direct() function.
                    Returns nothing.
      model       - model variable containing data needed to run simulation
      measurements - measurements structure. This value is also passed to the
                    'on_measurement' function (if specified).
      residual_fn - function with signature f(t, z, zdot, result, model)
                    where:
                        t - current time
                        z - current variables vector
                        zdot - vector of derivatives
                        results - vector of residuals, needs to be filled
                        model - supplied model variable
                    Returns nothing.
      initialize_zp0 - function with signature f(zp0, model). This allows to
                    provide the solver with better estimates of derivatives,
                    which in turn can speed-up the starting of the computation.
                    Function arguments are:
                        zp0 - vector to be filled with initial values
                              of derivative (or it's estimates)
                        model - argument supplied to simulate_direct() function.
                    Returns nothing.
      update_initial_condition - function with signature f(z0, model) or None.
                    The function is called on phase changing and allows to
                    adjust the re-starting condition for the solver.
                    Function arguments are:
                        z0 - vector of initial values of the next phase
                        model - supplied model variable
                    Returns nothing
      root_fn     - rooting function with signature f(t, z, zdot, result, model)
                    or None. If the function is present, the solver stops when
                    for any 'i' the value of result[i] is zero. Depending on the
                    'continue_on_root_found' argument the solver handles further
                    computation.
                    Function arguments:
                        t - current time
                        z - current variables vector
                        zdot - vector of derivatives
                        results - vector of values that needs to be filled;
                            the vector is of length nr_rootfns (see also
                            'nr_rootfns' argument). If some value
                            in results is zero, computation terminates.
                            (see also 'continue_on_root_found' argument)
                        model - supplied model variable
                    Returns nothing.
      nr_rootfns   - number of rooting functions (length of results vector, see
                    also the 'root_fn' argument). Value: interger.
      continue_on_root_found - function with signature f(model, t, z) or None.
                    Function can decide to return True or False depending on the
                    supplied arguments. If False is returned, computation is
                    aborted and the solver returns values found so far.
                    Function arguments:
                        model - supplied model variable
                        t - current time
                        z - current variables (values) vector
                    Returns True/False.
      algvars_idx  - vector or None. If vector, then specifies the indexes of
                    the vector 'z' which value is specified as algebraic
                    condition instead of differential one.
      on_phase_change - None or function with signature f(model, phase). This
                    function is called when phases are changed to allow adapt
                    the internal state of supplied data - model.
                    Function arguments:
                        model - supplied model variable
                        phase - new phase; phase is one of:
                            'a' - acceleration/centrifugation phase
                            'd' - deceleration phase
                            'g' - gravitation only
                    Returns nothing.
      on_measurement - None or a function with signature
                    f(t, z, model, measurements). The function is called when
                    a time of measurement is reached.
                    Function arguments:
                        t - current time
                        z - current variables vector
                        model - supplied model variable
                        measurements - the measurements object (see also the
                            'measurements' argument).
                    Returns nothing.
      truncate_results_on_stop - Type: boolean
                    If true, on root found or error returns only computed
                    data and values at times not reached will be truncated
                    (by default an array of size for all results is
                    preallocated). Otherwise an full array containing also
                    not-reached times data is returned.
    Return values:
      flag - Type: boolean
           if True, computation was sucessfully finished; otherwise
           a root was found/an error occured and value is set to False
      t  - returned times at which solver sucessfully computed results;
           if an error was encountered, the time at which it happened is
           stored in t[i] (for index i see below);
           if truncate_results_on_stop=False, remaining values are 0.0
      z  - returned sucessfully computed values; if an error occured,
           then z[i, :] (see index i below) contains the value at the
           time the error occurred
           if truncate_results_on_stop=False, remaining values are 0.0
      i  - index of last sucessful measurement (returned only if
           truncate_results_on_stop = False)
    """

    # Check supplied arguments
    if (not root_fn is None) and (nr_rootfns is None):
        raise Exception("Error: Function 'root_fn' was set in "
                        "simulate_inverse, but not the number of rooting "
                        "functions 'nr_rootfns'.")

    # Initialization
    measurements.reset_calc_measurements() # reset if previously stored values

    verbosity = model.verbosity

    measurements_times = measurements.get_times()

    measurements_nr = np.alen(measurements_times)

    if truncate_results_on_stop:
        t   = np.empty([measurements_nr, ], dtype=float)
        z   = np.empty([measurements_nr, model.z_size], dtype=float)
    else:
        t   = np.zeros([measurements_nr, ], dtype=float)
        z   = np.zeros([measurements_nr, model.z_size], dtype=float)

    t0 = 0.0

    model.set_value('omega_start', 0.0)
    model.set_value('t0', t0)
    model.set_value('phase', 'U') # initialize to dummy value

    z0 = np.empty([model.z_size, ], float)
    zp0 = np.zeros(z0.shape, float)

    initialize_z0(z0, model) # compute the initial state
    if not initialize_zp0 is None: initialize_zp0(zp0, z0, model)

    if measurements_times[0] == t0:
        t[0]    = t0
        z[0, :] = z0
        i = 1
        if not on_measurement is None:
            on_measurement(t0, z0, model, measurements)
    else:
        i = 0
    t_meas = measurements_times[i]

    solver = ida.IDA(residual_fn,
                     old_api=False,
                     compute_initcond=model.compute_initcond,
                     first_step_size=model.first_step_size,
                     atol=model.atol, rtol=model.rtol,
                     max_step_size = model.max_step_size,
                     max_steps = model.max_steps,
                     algebraic_vars_idx=algvars_idx,
                     exclude_algvar_from_error=exclude_algvar_from_error,
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
            print ('adding tstop at:', t_end, 'Phase:', phase)

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
                solution = solver.init_step(t0, z0, zp0)
                if solution.errors.t:
                    print ('ERROR:', solution.message)
                    if truncate_results_on_stop:
                        return (False, t[:i+1], z[:i+1, :], i-1)
                    else:
                        return (False, t, z, i-1)
                #t_init = solution.values.t

                solver_initialized = True

            if not on_phase_change is None:
                if previous_phase is None: previous_phase = phase
                if previous_phase != phase: on_phase_change(model, phase)

            perc_done=0;
            perc_next=1;
            while True:
                #print ( ' Step start', t_meas )
                #(flag, t_out) = solver.step(t_meas, z[i, :], zp0)
                solution = solver.step(t_meas, z[i, :], zp0)

                if solution.errors.t:     # error occured
                    if verbosity > 1:
                        print('Error occured during computation. Solver failed '
                              'at time\nt_err=', solution.errors.t,
                              '\nExpected value of t:', t)
                    if verbosity > 2:
                        print('Values at this time:\nz_err=', z[i, :])

                    if truncate_results_on_stop:
                        return (False, t[:i+1], z[:i+1, :], i)
                    else:
                        t[i] = solution.errors.t
                        return (False, t, z, i)
                elif solution.roots.t: # root found
                    t[i] = solution.values.t
                    if ((not continue_on_root_found is None)
                         and (continue_on_root_found(model, t[i], z[i, :]))):
                        if verbosity > 1:
                             print('Root found: continuing computation')
                        if truncate_results_on_stop:
                            return (False, t[:i+1], z[:i+1, :], i)
                        else:
                            return (False, t, z, i)
                    else:
                        if verbosity > 1:
                             print('Root found: aborted further computation.')
                        if not on_measurement is None:
                            on_measurement(t[i], z[i, :], model,
                                           measurements)
                        if truncate_results_on_stop:
                            return (True, t[:i+1], z[:i+1, :], i)
                        else:
                            return (True, t, z, i)

                else: # success or t_stop reached
                    t[i] = solution.values.t
                    if (solution.flag == StatusEnumIDA.SUCCESS) or (t[i] == t_meas):
                        if not on_measurement is None:
                            on_measurement(t[i], z[i, :], model,
                                           measurements)

                        i += 1
                        if i < measurements_nr:
                            t_meas = measurements_times[i]
                        else:
                            # previous measurement was the last measurement we
                            # have taken so we can abort further computation
                            return (True, t, z, i-1)

                    # Assuming t_out is t_end+numerical_error
                    if (solution.flag == StatusEnumIDA.TSTOP_RETURN) or (t_end == t[i]):
                        break

                perc_done = (t[i-1]-model.t0) / (t_end-model.t0) * 100
                if perc_done >= perc_next:
                    if show_progress:
                        print ('Finished %03d percent of computation' % int(perc_done))
                    perc_next = int(perc_done+1)

            t0 = t_end

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
            if not np.all(z0 == z[i-1, :]):
                solver_initialized = False
        else:
            z0 = z[i-1, :]

        if model.always_restart_solver:
            #makes sense only if there is an init cond to compute
            if model.compute_initcond == None:
                print ("Error: Asking to restart solver, but no computation of init condition")
                sys.exit(0)
            solver_initialized = False

    return (True, t, z, i-1)

def print_params(params):
    """ Print supplied parameters of type dict. """

    print()
    for (name, value) in params.items():
        if name == 'ks':
            if np.iterable(value):
                print('{0:5}: {1!s}'.format('Ks', value))
            else:
                print('{:5}: {: .8g}'.format('Ks', value))
        else:
            if np.iterable(value):
                print('{0:5}: {1!s}'.format(name, value))
            else:
                print('{:5}: {: .8g}'.format(name, value))

def untransform_dict(names, values, lengths, result, untransform):
    inilen = 0
    if untransform:
        for (name, length) in zip(names, lengths):
            result[name] = untransform[name](values[inilen:inilen+length])
            inilen += length
    else:
        for (name, length) in zip(names, lengths):
            result[name] = values[inilen:inilen+length]
            inilen += length

def penalize(parameters, lbounds, ubounds, when='out_of_bounds'):
    """
      Compute penalization coeficient.
      If 'when' is set to 'out_of_bounds', the coef is found a if boundaries
      are crossed (if not, 0 is returned); otherwise a coef based on the
      distance of parameters from boundary is computed.
    """
    max_penalization = 1e50

    penalization = 0.0

    if when == 'out_of_bounds':
        tolerance = 0.15
        #a_max = 0.01
        for (name, values) in parameters.items():
            if not np.iterable(values):
                values = [values]
            for value in values:
                if lbounds[name] > value:
                    a = np.abs(value - lbounds[name])
                elif ubounds[name] < value:
                    a = np.abs(value - ubounds[name])
                else:
                    continue

                penalty = (1/tolerance) + (10*a)

                penalization += penalty # min(penalty, max_penalization)
                print('penalty', penalty, penalization)
    else:
        #failed computation
        for (name, values) in parameters.items():
            if not np.iterable(values):
                values = [values]
            for value in values:
                a = min(np.exp(value - lbounds[name]), max_penalization)
                b = min(np.exp(value - ubounds[name]), max_penalization)

                penalization += 1e3*min(10 * (a + 1/a) + 10 * (b + 1/b),
                                    max_penalization)

    return penalization

def penalize_cond(parameters, conditions):

    #max_penalization = 1e5
    tolerance = 0.15

    penalization = 0.0
    for (name, values) in parameters.items():
        if not np.iterable(values):
            values = [values]
        for cond in conditions[name]:
            if cond == '':
                continue
            if cond.lower() in ['ascending', 'asc']:
                #parameters must be ascending
                prevval = values[0]
                for value in values[1:]:
                    if value < prevval:
                        a = prevval - value
                    else:
                        prevval = value
                        continue
                    penalty = 1/(tolerance) + 50*a
                    penalization += penalty # min(penalty, max_penalization)
            elif cond.lower() in ['positive', 'pos']:
                for value in values:
                    if value <= 0:
                        penalty = 1/(tolerance) + 50*(-value)
                        penalization += penalty # min(penalty, max_penalization)
            else:
                raise Exception('Condition on parameters not recoginized: %s' % cond)
    return penalization

def simulate_inverse(direct_fn, model, measurements, optimfn='leastsq'):
    """
    As long as refinements in approach are possible, restart inverse method
    """
    invrun =0;

    while True:
        (optim_params, cov) = run_inverse(direct_fn, model, measurements, optimfn)
        invrun += 1
        # make sure we work with optim parameters
        print ('setting optim params for final run', optim_params)
        model.set_parameters(optim_params)

        # reset if previously stored values of measurements (it's re-set also
        # inside the simulate_direct() but user may supply his own direct
        # function, so we need to be sure measurements were re-setted
        measurements.reset_calc_measurements()

        #output inbetween data of this refinement run.
        show_inbetween_result(model.experiment_info, model,
                              optim_params, cov, ext='_%03d' % invrun)

        print ("Doing refine now")

        if hasattr(model, "refine"):
            success = model.refine()
            print ("Doing refine .. ", success)
            if not success:
                print ("Stopping refinement, current par")
                print (model.get_parameters(optim_params))
                break;
            else:
                #refine was success, we update the params
                print ('Refinement success, current par')
                print (model.get_parameters(optim_params))
        else:
            break;

    return (optim_params, cov)

def run_inverse(direct_fn, model, measurements, optimfn='leastsq'):
    available_solvers = ['leastsq', 'fmin', 'fmin_powell', 'fmin_cg',
                         'fmin_bfgs', 'raster', 'fmin_slsqp']
    if not optimfn in available_solvers:
        print("Unknown inverse method solver for 'optimfn': ", optimfn)
        print("Available solvers are: ", available_solvers)
        exit(1)

    untransform  = model._untransform
    lbounds      = model._lbounds
    ubounds      = model._ubounds
    conditions   = model._conditions
    optim_params = {}
    global last_good_res
    last_good_res = None

    global ITERATION # Python 2.7 hack (no support for nonlocal variables)
    ITERATION = 0

    def optimfn_wrapper(optimargs):
        global ITERATION # Python 2.7 hack (no support for nonlocal variables)
        global last_good_res

        print(15 * '*', ' Iteration: {:4d}'.format(ITERATION), ' ', 15 * '*')
        ITERATION += 1

        print ( "True arg", optimargs)
        untransform_dict(optim_names, optimargs, optim_par_length,
                         optim_params, untransform)

        if model.verbosity > 0: print_params(optim_params)

        penalization = penalize(optim_params, lbounds, ubounds,
                                when='out_of_bounds')
        penalization += penalize_cond(optim_params, conditions)
#        if penalization > 0:
#            print ('not ascending par', penalization)

        if penalization > 0.0:
            if penalization < 1e3: penalization += 1e3
            if model.verbosity > 0:
                print('Optimized arguments are out of bounds or not satisfying conditions ... Penalizing by ',
                      penalization)

            #error = measurements.get_penalized(penalization,
            #                scalar=(optimfn not in ['leastsq']))
            if optimfn in ['leastsq']:
                error = np.abs(last_good_res) + penalization/len(last_good_res)
            else:
                error = np.abs(last_good_res) + penalization
            #if model.verbosity > 0:
            #    print('Penalized error is', error)
            return error

        print ( "Used arg", optim_params)
        model.set_parameters(optim_params)

        # reset if previously stored values of measurements (it's re-set also
        # inside the simulate_direct() but user may supply his own direct
        # function, so we need to be sure measurements were re-setted
        measurements.reset_calc_measurements()

        flag = direct_fn(model, measurements)

        if flag:
            # direct computation went O.K.
            if model.verbosity > 0:
                measurements.display_error()

            error = measurements.get_error() # computed - measured

            total_LSQ_error = np.sum(np.power(error, 2))
            print('\nTotal LSQ error:', total_LSQ_error)

            if optimfn in ['leastsq']:
                result = error
            else:
                result = total_LSQ_error

            last_good_res = result

        else:
            # something is wrong, so penalize
            penalization = \
              min(penalize(optim_params, lbounds, ubounds, when='always'),
                  1e10)
            if model.verbosity > 1:
                print('Direct problem did not converge for given optimization '
                      'parameters... Penalizing by ', penalization)

            result = measurements.get_penalized(penalization,
                            scalar=(optimfn not in ['leastsq']))

        return result

    def optimineq_wrapper(optimargs):
        #inequality conditions

        untransform_dict(optim_names, optimargs, optim_par_length,
                         optim_params, untransform)
        ineqs = np.empty((0,), float)
        for (name, values) in optim_params.items():
            if not np.iterable(values):
                values = [values]
            for cond in conditions[name]:
                if cond == '':
                    continue
                if cond.lower() in ['ascending', 'asc']:
                    #parameters must be ascending
                    ineqs = np.append(ineqs, values[1:]-values[:-1])
                else:
                    raise Exception('Condition on parameters not recoginized: %s' % cond)
        if len(ineqs) == 0:
            raise Exception('Use inequality solver only if there are are inequality'
                            'conditions! None found.')
        return ineqs

    # Further it is assumed that we are using OrderedDict so that the
    # order of values is stable
    from collections import OrderedDict
    assert type(model.init_values) is OrderedDict, \
        'The type of ''model.init_values'' is not OrderedDict.'

    optim_names  = model.init_values.keys()   # it's an ordered dict
    optim_par_length = np.ones((len(optim_names), ), dtype=int)

    # Determine lengths of each parameter (in case value is an array)
    for (ind, val) in enumerate(model.init_values.values()):
        if np.iterable(val):
            optim_par_length[ind] = len(val)

    transform = model._transform

    if bool(transform is None) == bool(untransform is None):
        if not transform is None:
            from .saturation_curve import default_transformation

            for name in optim_names:
                # Add missing [un]transform functions
                if bool(name in transform) == bool(name in untransform):
                    if not name in transform:
                        (transform[name], untransform[name]) = \
                          default_transformation(-np.inf, np.inf)
                else:
                    raise Exception('Name ''', name, ' '' is specified in one of '
                        '[un]transform but not in the other:\n'
                        'Transform:   ', transform,
                        '\nUntransform: ', untransform)
    else:
        raise Exception('One of transform/untransform is ''None'' while the other not.')

    init_values = np.empty((np.sum(optim_par_length), ), dtype=float)
    iv_ind = 0

    for (ind, name) in enumerate(optim_names):
        # Update init_values
        iv_ind_next = iv_ind + optim_par_length[ind]
        if transform and name in transform:
            init_values[iv_ind:iv_ind_next] = \
            transform[name](np.asarray(model.init_values[name]))
        else:
            init_values[iv_ind:iv_ind_next] = np.asarray(model.init_values[name])
        iv_ind = iv_ind_next

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
    for run in range(1 + (model.transform_params and model.untransformed_cov)):
        # if untransformed covariance is desired but we computed with
        # transformed parameters, we run the computation again, starting
        # in optimal value of previous optimization
        if run == 1:
            print('Running optimization to obtain the covariance of '
                  'untransformed parameters...')
            untransform = None
            init_values = list(flatten([optim_params[name] for name in optim_names]))

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
        elif optimfn == 'fmin_slsqp':
            (opt_params, fopt,fcalls, warnflag, msg) = \
                optimize(optimfn_wrapper, init_values,
                         f_eqcons=optimfn_wrapper,
                         f_ieqcons=optimineq_wrapper,
                         iter=model.max_inv_iter or 100,
                         acc=model.xtol, full_output=True,
                         epsilon=model.epsfcn, disp=model.disp_inv_conv
                         )
        print ("Final param", opt_params)
        untransform_dict(optim_names, opt_params, optim_par_length,
                         optim_params, untransform)
        print ("unstransf final param", optim_params)

    # run the direct once more to obtain correct covariance
    print ('rerunning for opt error')
    opt_error = optimfn_wrapper(opt_params)
    print ('finished, error:', opt_error)

    if cov is None:
        print('Warning: singular matrix for covariance  encountered '
              '(indicates very flat curvature in some direction).')
    else:
        s_sq = (np.sum(np.power(opt_error, 2))
                / (np.alen(opt_error) - np.alen(opt_params)))
        cov *= s_sq

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
    out += ' |\n'
    for (name, value) in results:
        if not value is None:
            out += ' |{:8d}'.format(value)
    print(out, '|')
    print("Extra info:", opt_params, cov, infodic, msg, ier)

    return (optim_params, cov)

#def compute_raster(model):
#        lbounds = xdata.inv_lbounds
#        ubounds = xdata.inv_ubounds
#        raster_size = xdata.raster_grid_size
#        nrlevels = len(lbounds)
#        if np.isscalar(raster_size):
#            raster_size = [raster_size] * nrlevels
#
#        optimargs = [0]* nrlevels
#        grid = []
#        output = []
#
#        from copy import copy
#        for (lb, ub, rsize) in zip(lbounds, ubounds, raster_size):
#            grid.append(np.linspace(float(lb), float(ub), float(rsize)))
#
#        def loopsingle(level):
#            for (ind, val) in enumerate(grid[level]):
#                optimargs[level] = val
#
#                if level == nrlevels - 1:
#                    output.append((fmin_wrapper_fn(optimargs, xdata),
#                                  copy(optimargs)))
#                else:
#                    loopsingle(level+1)
#
#        loopsingle(0)
#
#        print(output)
#        from  pickle import dump
#        dump(output, '/home/archetyp/Desktop/dump.bin')
