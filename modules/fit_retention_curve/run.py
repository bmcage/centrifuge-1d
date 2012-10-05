import numpy as np
from scipy.optimize import leastsq
from modules.shared.vangenuchten import h2u
from modules.shared.functions import show_results

def solve(model):
    theta_measured = model.measurements['theta'][0]

    def direct_wrapper(optim_args):
        if model.verbosity > 1:
            for (name, value) in zip(params_names, optim_args):
                print('{}: {: >10.7f}     '.format(name, value), end='')
            print()

        for (name, value) in zip(params_names, optim_args):
            setattr(model, name, value)
            if name == 'n':
                setattr(model, 'm', 1. - 1./value)

        u = h2u(model.h, model.n, model.m, model.gamma)

        (theta_s, theta_r) = (model.theta_s, model.theta_r)
        theta = theta_r + (theta_s - theta_r) * u

        return np.power(theta - theta_measured, 2)

    #    backup_params = model.get_parameters(model.inv_init_params.keys())

    # [p] = Pa = kg/m/s^2 = 10 * g/cm/s^2 -\
    # [h] = cm                            - \
    # => h [cm] =(10*p)/g/rho with [g]=cm/s^2, [rho]=g/cm^3
    p = np.asarray(model.p, dtype=float)
    model.h = -10.*p / model.rho /model.g

    (params_names, init_values) = ([], [])
    for (name, value) in model.inv_init_params.items():
        params_names.append(name)
        init_values.append(value)

    (opt_values, cov, infodic, msg, ier) = \
          leastsq(direct_wrapper, init_values,
                   epsfcn=model.epsfcn, factor=model.factor,
                   xtol=model.xtol, ftol=model.ftol,
                   full_output=True)

    optim_params = {key: value
                    for (key, value) in zip(params_names, opt_values)}

    return (optim_params, cov)

P_DISP = np.arange(0, 10000000, 100)

def extract_data(model):
    global P_DISP
    (n, gamma, theta_s, theta_r) = (model.n, model.gamma,
                                    model.theta_s, model.theta_r)
    theta = theta_r + ((theta_s - theta_r)
                       * h2u(-10.*P_DISP/model.rho/model.g, n, 1.-1./n, gamma))
    theta_in_measured_points = theta_r + ((theta_s - theta_r)
                                          * h2u(model.h, n, 1.-1./n, gamma))

    extracted_data = {'theta': (theta, P_DISP, theta_in_measured_points)}

    return (True, extracted_data)

def run(model):
    (inv_params, cov) = solve(model)

    # DISPLAY RESULTS:
    if inv_params:
        model.set_parameters(inv_params)
        model_verbosity = model.verbosity # backup verbosity
        model.verbosity = 0

        show_results(extract_data, model, inv_params=inv_params, cov=cov)

        model.verbosity = model_verbosity # restore verbosity

    return inv_params
