import scikits.odes.sundials.ida as ida
from numpy import (zeros, concatenate, all, sum, power, isscalar, linspace,
                   asarray, cumsum, min)

simulation_err_str = ('%s simulation: Calculation did not reach '
                      'the expected time. An error occured.'
                      'Reached time: % 10.6f\nExpected time: % 10.6f')

class DirectSimulator:
    def __init__(self, model, residual_fn, algvars_idx = None, root_fn = None,
                 nr_rootfns=None):

        self.model = model
        self.residual_fn = residual_fn

        self.root_fn = root_fn
        if (not root_fn is None) and (nr_rootfns is None):
            raise Exception("Error: Function 'root_fn' was set in DirectSimulator, "
                            "but the number of root 'nr_rootfns' not.")
        else:
            self.nr_rootfns = nr_rootfns
        self.algvars_idx = algvars_idx

        self.z_previous  = None
        self.zp_previous = None
        self.solver = None
        self.t0 = 0.0

    def run(self, duration, fh_duration, z0, z_retn):
        if self.z_previous is None:
            model = self.model
            self.solver = ida.IDA(self.residual_fn,
                                  compute_initcond='yp0',
                                  first_step_size=model.first_step_size,
                                  atol=model.atol, rtol=model.rtol,
                                  max_step_size = model.max_step_size,
                                  max_steps = model.max_steps,
                                  algebraic_vars_idx=self.algvars_idx,
                                  rootfn=self.root_fn,
                                  nr_rootfns=self.nr_rootfns,
                                  tcrit=duration,
                                  #linsolver='band', uband=1, lband=1,
                                  user_data=model)

            self.z_previous = zeros(z0.shape, float)

        verbosity = self.model.verbosity

        if self.model.always_restart_solver or (not all(self.z_previous == z0)):
            zp0 = zeros(z0.shape, float)
            self.solver.init_step(self.t0, z0, zp0)
        if duration > 0.0:
            t = self.t0 + duration

            if verbosity > 1:
                if self.t0 == 0.0:
                    print(27*' ', end='')
                print(' ... Starting time: %9.1f  ' % self.t0,
                      'Duration: %9.1f  ' % duration,
                      'Expected end time: %9.1f' % t)

            self.solver.set_options(tcrit=t)
            flag, t_out = self.solver.step(t, z_retn)

            if t_out < t:
                print('Error occured during computation. Solver returned with '
                      'values:\nt_err=', t_out, '\nz_err=', z_retn,
                      '\nExpected value of t:', t)
                return (False, t_out)

            if verbosity > 1:
                print('Computed end time: %8.1f' % t_out, end='')
            self.t0 = t_out
            self.z_previous[:] = z_retn

        if fh_duration > 0.0:
            t = self.t0 + fh_duration
            self.solver.set_options(tcrit=t)

            acceleration = model.include_acceleration
            r0           = model.r0
            omega        = model.omega
            duration     = model.duration

            model.include_acceleration = False
            model.r0    = model.r0_fall
            model.omega = model.omega_fall
            model.duration = fh_duration

            flag, t_out = self.solver.step(t, z_retn)

            if flag < 0:
                print('Error occured during computation. Solver returned with '
                      'values:\nt_err=', t_out, '\nz_err=', z_retn,
                      '\nExpected value of t:', t)
                return (False, t_out)

            self.t0 = t_out
            self.z_previous[:] = z_retn

            model.include_acceleration = acceleration
            model.r0    = r0
            model.omega = omega
            model.duration = duration

        return (True, t_out)


def simulate_direct(model, residual_fn, z0, algvars_idx = None, root_fn = None,
                    nr_rootfns=None, first_step_size=1e-30):

    def run_simulation(duration, z0):
        solver = ida.IDA(residual_fn,
                         compute_initcond='yp0',
                         first_step_size=first_step_size,
                         atol=model.atol, rtol=model.rtol,
                         max_step_size = model.max_step_size,
                         max_steps = model.max_steps,
                         algebraic_vars_idx=algvars_idx,
                         rootfn=root_fn, nr_rootfns=nr_rootfns,
                         tcrit=duration,
                         #linsolver='band', uband=1, lband=1,
                         user_data=model)

        zp0 = zeros(z0.shape, float)

        flag, t, z, zp, t_err, z_err, zp_err = solver.solve([0.0, duration], z0, zp0)

        if flag == 0: # Success
            return flag, t[1], z[1]
        elif flag == 1: # TCrit reached,also success
            return flag, t_err, z_err
        elif flag == 2: # Root found, considered to be an success
            return flag, t_err, z_err
        else:
            print('Error occured during computation. Solver returned with '
                  'values:\nt_err=', t_err, '\nz_err=', z_err, '\nzp_err=',
                  zp_err)
            exit(1)

    t_retn_out = 0.0

    if model.duration > 0.:
        duration = model.duration
        flag, t_out, z_out = run_simulation(duration, z0)

        if t_out < duration:
            print(simulation_err_str % ('Centrifugation', t_out, duration))
            return (False, t_out, z_out)
        z0 = z_out
        t_retn_out = t_out

    if model.fh_duration > 0.:
        acceleration = model.include_acceleration

        model.include_acceleration = False
        model.r0    = model.r0_fall
        model.omega = model.omega_fall

        duration = model.fh_duration
        flag, t_fh_out, z_out = run_simulation(duration, z0)

        model.include_acceleration = acceleration

        t_retn_out = t_retn_out + t_fh_out

        if t_fh_out < duration:
            print(simulation_err_str % ('Falling head test',
                                        t_fh_out, duration))
            return (False, t_retn_out, z_out)

    return (True, t_retn_out, z_out)

def simulate_inverse(times, direct_fn, model, init_parameters,
                     wl_in_meas=None, wl_out_meas=None, gc_meas=None,
                     rm_meas=None,
                     optimfn='leastsq'):

    from modules.shared.functions import determine_scaling_factor
    from modules.shared.show import disp_inv_results
    from numpy import empty, log, exp, alen

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
            for name in optimized_parameters:
                print('%5s: %f' % (name, getattr(model, name)))

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
            for (param, value) in zip(optimized_parameters, optim_args):
                a = exp(value - lbounds[param])
                b = exp(value - ubounds[param])

                penalization = (penalization + 10 * (a + 1/a) + 10 * (b + 1/b))

        return penalization

    calc_wl_in  = bool(wl_in_meas)
    calc_wl_out = bool(wl_out_meas)
    calc_gc     = bool(gc_meas)
    calc_rm     = bool(rm_meas)

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

    def optimfn_wrapper(optimargs):
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
                                 display_graphs=False)


            computation = concatenate((wl_in_C, wl_out_C, gc_C, rm_C))

            if optimfn == 'leastsq':
                return (computation - measurements)
            else:
                return sum(power(computation - measurements, 2))

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
    (flag, t, wl_in, wl_out, gc, rm) = direct_fn(model)
    disp_inv_results(model, t, inv_params=optim_params, cov=cov,
                     wl_in_inv=wl_in, wl_out_inv=wl_out, gc1_inv=gc, rm1_inv=rm)

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
