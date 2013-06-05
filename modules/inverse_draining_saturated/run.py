from __future__ import division
from __future__ import print_function
from modules.direct_draining_saturated.run import solve as solve_direct
from modules.shared.show import show_results
from modules.shared.solver import simulate_inverse

def solve(model, measurements):

    (inv_params, cov) = \
      simulate_inverse(solve_direct, model, measurements, optimfn=model.optimfn)

    return (inv_params, cov)

def run(model):
    (inv_params, cov) = solve(model, model.measurements)

    # DISPLAY RESULTS:
    if inv_params:
        model.set_parameters(inv_params)
        # run once again the direct problem with optimal parameters
        model_verbosity = model.verbosity # backup verbosity
        model.verbosity = 0

        show_results(model.experiment_info, model=model,
                     inv_params=inv_params, cov=cov)

        model.verbosity = model_verbosity # restore verbosity

    return inv_params

def dry_run(model):
    import numpy as np

    inv_params = {'n': 2.81, 'gamma': -0.0189, 'ks': 1e-5}
    cov = np.asarray([[2.0, 3.0, 4.0], [5.0, 6.0, 7.0], [1.0, 3.0, 6.0]])

    # DISPLAY RESULTS:
    if inv_params:
        model.set_parameters(inv_params)
        # run once again the direct problem with optimal parameters
        model_verbosity = model.verbosity # backup verbosity
        model.verbosity = 0

        show_results(model.experiment_info, model=model,
                     inv_params=inv_params, cov=cov)

        model.verbosity = model_verbosity # restore verbosity

    return inv_params
