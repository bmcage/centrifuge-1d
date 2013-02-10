from __future__ import division
from __future__ import print_function
from modules.direct_draining_saturated.run import solve as solve_dds, \
     extract_data
from modules.shared.functions import show_results
from modules.shared.solver import simulate_inverse

from modules.shared.vangenuchten import h2u
from numpy import alen

def solve(model, measurements):

    #calc_p = model.get_parameters(('calc_gc', 'calc_rm', 'calc_wm', 'calc_f_mt',
    #                               'calc_f_mo', 'calc_cf_mo')) # backup
    (inv_params, cov) = \
      simulate_inverse(solve_dds, model, model.inv_init_params,
                       measurements, optimfn=model.optimfn)

    #model.set_parameters(calc_p) # ...and restore values

    return (inv_params, cov)

def run(model):
    (inv_params, cov) = solve(model, model.measurements)

    # DISPLAY RESULTS:
    if inv_params:
        model.set_parameters(inv_params)
        # run once again the direct problem with optimal parameters
        model.calc_wm = True
        model_verbosity = model.verbosity # backup verbosity
        model.verbosity = 0

        show_results(model.experiment_info, model=model,
                     inv_params=inv_params, cov=cov)

        model.verbosity = model_verbosity # restore verbosity

    return inv_params
