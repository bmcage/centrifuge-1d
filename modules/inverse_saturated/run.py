from __future__ import division

from modules.direct_saturated.run import solve as direct_solve, extract_data
from modules.shared.solver import simulate_inverse

def solve(model):

    def ip_direct_saturated_heights(model, measurements_names):

        (flag, t, z, measurements) = direct_solve(model)

        result = [flag, t]

        for name in measurements_names:
            result.append(measurements[name][1:]

        return result

    inv_params = simulate_inverse(ip_direct_saturated_heights, model,
                                  model.inv_init_params, optimfn=model.optimfn)
    return inv_params

def run(model):
    (inv_params, cov) = solve(model)

    if inv_params:
        from modules.shared.functions import show_results

        model.set_parameters(inv_params)
        # run once again the direct problem with optimal parameters
        model.calc_wm = True
        model_verbosity = model.verbosity # backup verbosity
        model.verbosity = 0

        show_results(model.experiment_info, model=model,
                     inv_params=inv_params, cov=cov)

        model.verbosity = model_verbosity # restore verbosity
