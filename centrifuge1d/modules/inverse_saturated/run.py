from __future__ import division

from ..direct_saturated.run import solve as solve_direct
from ..shared.solver import simulate_inverse

def solve(model, measurements):

    (inv_params, cov) = \
      simulate_inverse(solve_direct, model, measurements, optimfn=model.optimfn)

    return (inv_params, cov)

def run(model):
    (inv_params, cov) = solve(model, model.measurements)

    if inv_params:
        from ..shared.show import show_results

        model.set_parameters(inv_params)
        # run once again the direct problem with optimal parameters
        model.calc_wm = True
        model_verbosity = model.verbosity # backup verbosity
        model.verbosity = 0

        show_results(model.experiment_info, model=model,
                     inv_params=inv_params, cov=cov)

        model.verbosity = model_verbosity # restore verbosity
