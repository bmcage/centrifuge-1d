from __future__ import division

import numpy as np
from modules.shared.show import show_results
from modules.shared.solver import simulate_inverse

def solve_direct(model, measurements):
    # This is not needed when calling simulate_direct() as it resets the
    # values automatically
    measurements.reset_calc_measurements() # reset if previously stored values

    measurements.store_calc_theta(model.h, model.SC, model.theta_s,
                                  model.theta_r, model.rho, model.g)

    return True

def solve(model, measurements):
    # [p] = Pa = kg/m/s^2 = 10 * g/cm/s^2 -\
    # [h] = cm                            - \
    # => h [cm] =(10*p)/g/rho with [g]=cm/s^2, [rho]=g/cm^3
    p = np.asarray(model.p, dtype=float)
    model.h = -10.*p / model.rho /model.g

    (inv_params, cov) = \
      simulate_inverse(solve_direct, model, measurements, optimfn=model.optimfn)

    return (inv_params, cov)

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
