import scikits.odes.sundials.ida as ida
from numpy import zeros
from scipy.optimize import curve_fit

simulation_err_str = ('%s simulation: Calculation did not reach '
                      'the expected time. An error occured.'
                      'Reached time: % 10.6f\nExpected time: % 10.6f')

def simulate(model, residual_fn, z0):

    def run_simulation(duration, z0):
        solver = ida.IDA(residual_fn,
                         compute_initcond='yp0',
                         first_step_size=1e-18,
                         atol=model.atol, rtol=model.rtol,
                         max_step_size = 840.,
                         max_steps = 8000,
                         #algebraic_vars_idx=[4],
                         linsolver='band', uband=1, lband=1,
                         user_data=model)

        zp0 = zeros(z0.shape, float)

        flag, t, z = solver.solve([0.0, duration], z0, zp0)[:3]

        return flag, t[1], z[1]


    if model.duration > 0.:
        duration = model.duration
        flag, t_out, z_out = run_simulation(duration, z0)

        if t_out < duration:
            print(simulation_err_str % ('Centrifugation', t_out, duration))
            return (False, t_out, z_out)
        z0 = z_out

    elif model.fh_duration > 0.:
        acceleration = model.include_acceleration

        model.include_acceleration = False
        model.r0    = model.r0_fall
        model.omega = model.omega_fall

        duration = model.fh_duration
        flag, _t_out, z_out[:]= run_simulation(duration, z0)

        model.include_acceleration = acceleration

        if t_out < duration:
            print(simulation_err_str % ('Falling head test', t_out, duration))
            return (False, t_out, z_out)

    return (True, t_out, z_out)

def simulate_inverse(lsq_direct_fn, xdata, ydata, init_params):
    return curve_fit(lsq_direct_fn, xdata, ydata, p0 = init_params)
