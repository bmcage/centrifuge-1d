from __future__ import division

from modules.direct_saturated.run import solve as solve_ds, extract_data
from modules.shared.solver import simulate_inverse

def solve(model, measurements):

    inv_params = simulate_inverse(solve_ds, model, model.inv_init_params,
                                  measurements, optimfn=model.optimfn)

    return inv_params

def run(model):
    (inv_params, cov) = solve(model, model.measurements)

    if inv_params:
        from modules.shared.show import show_results

        model.set_parameters(inv_params)
        # run once again the direct problem with optimal parameters
        model.calc_wm = True
        model_verbosity = model.verbosity # backup verbosity
        model.verbosity = 0

        show_results(model.experiment_info, model=model,
                     inv_params=inv_params, cov=cov)

        model.verbosity = model_verbosity # restore verbosity
