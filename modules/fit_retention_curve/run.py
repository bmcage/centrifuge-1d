from __future__ import division

import numpy as np
from scipy.optimize import leastsq
from modules.shared.vangenuchten import h2u, retention_curve
from modules.shared.functions import show_results

def solve(model):
     # list of measured data - we have only one measurement, so first item
    theta_measured = model.measurements.get_values()[0]

    def direct_wrapper(optim_args):
        out = ''
        if model.verbosity > 1:
            for (name, value) in zip(params_names, optim_args):
                out += '{}: {: >10.7f}     '.format(name, value)
            print(out)

        for (name, value) in zip(params_names, optim_args):
            setattr(model, name, value)
            if name == 'n':
                setattr(model, 'm', 1. - 1./value)

        u = h2u(model.h, model.n, model.m, model.gamma)

        (theta_s, theta_r) = (model.theta_s, model.theta_r)
        theta = theta_r + (theta_s - theta_r) * u

        return theta - theta_measured

    #    backup_params = model.get_parameters(model.inv_init_params.keys())

    # [p] = Pa = kg/m/s^2 = 10 * g/cm/s^2 -\
    # [h] = cm                            - \
    # => h [cm] =(10*p)/g/rho with [g]=cm/s^2, [rho]=g/cm^3
    p = np.asarray(model.p, dtype=float)
    model.h = -10.*p / model.rho /model.g

    (params_names, init_values) = ([], [])
    for (name, value) in model.inv_init_params.items():
        params_names.append(name)
        init_values.append(value[0])

    (opt_values, cov, infodic, msg, ier) = \
          leastsq(direct_wrapper, init_values,
                   epsfcn=model.epsfcn, factor=model.factor,
                   xtol=model.xtol, ftol=model.ftol,
                   full_output=True)

    optim_params = {key: value
                    for (key, value) in zip(params_names, opt_values)}
    s_sq = np.power(direct_wrapper(opt_values),2).sum() / (len(model.h) - len(params_names))

    return (optim_params, cov * s_sq)

def extract_data(model):
    (n, gamma, theta_s, theta_r) = (model.n, model.gamma,
                                    model.theta_s, model.theta_r)

    (p, theta) = retention_curve(n, gamma, theta_s, model.rho,
                                 model.g, theta_r=theta_r)

    p_meas = np.asarray(model.p, dtype=float)
    (p_meas, theta_in_measured_points) = \
      retention_curve(n, gamma, theta_s, model.rho, model.g,
                      theta_r=theta_r, p=p_meas)

    extracted_data = {'theta': (p, theta, theta_in_measured_points)}

    return (True, extracted_data)

def run(model):
    (inv_params, cov) = solve(model)

    # DISPLAY RESULTS:
    if inv_params:
        model.set_parameters(inv_params)
        model_verbosity = model.verbosity # backup verbosity
        model.verbosity = 0

        show_results(model.experiment_info, model=model,
                     inv_params=inv_params, cov=cov)

        model.verbosity = model_verbosity # restore verbosity

    return inv_params
