import scikits.odes.sundials.ida as ida
from numpy import zeros, concatenate
from scipy.optimize import curve_fit

simulation_err_str = ('%s simulation: Calculation did not reach '
                      'the expected time. An error occured.'
                      'Reached time: % 10.6f\nExpected time: % 10.6f')

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

def simulate_inverse(lsq_direct_fn, xdata, ydata, init_params):

    def lsq_wrapper_fn(xdata, *optim_args):
        direct_results = lsq_direct_fn(xdata, optim_args)
        #print('direct results:', direct_results)
        return concatenate(direct_results)

    return curve_fit(lsq_wrapper_fn, xdata, ydata, p0 = init_params)
