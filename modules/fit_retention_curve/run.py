from __future__ import division

import numpy as np
from scipy.optimize import leastsq
from modules.shared.vangenuchten import h2u, retention_curve
from modules.shared.show import show_results
from modules.shared.solver import print_params

def solve(model, measurements):
     # list of measured data - we have only one measurement, so first item
    theta_measured = measurements.get_values()[0]

    def direct_wrapper(optim_args):
        params = dict(zip(params_names, optim_args))
        model.set_parameters(params)
        if model.verbosity > 1:
            print_params(params)

        (p, theta) = retention_curve(model.n, model.gamma, model.theta_s,
                                     model.rho, model.g, model.theta_r,
                                     p=None, h=model.h)

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

    optim_params = dict(zip(params_names, opt_values))
    s_sq = (sum(np.power(direct_wrapper(opt_values), 2))
            / (len(model.h) - len(params_names)))

    return (optim_params, cov * s_sq)

def run(model):
    (inv_params, cov) = solve(model, model.measurements)

    # DISPLAY RESULTS:
    if inv_params:
        model.set_parameters(inv_params)
        model_verbosity = model.verbosity # backup verbosity
        model.verbosity = 0

        show_results(model.experiment_info, model=model,
                     inv_params=inv_params, cov=cov)

        model.verbosity = model_verbosity # restore verbosity

    return inv_params
