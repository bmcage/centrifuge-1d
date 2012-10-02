from modules.direct_draining_saturated.run import solve as solve_direct, \
     extract_data
from modules.shared.functions import measurements_time, show_results
from modules.shared.solver import simulate_inverse
from modules.shared.vangenuchten import h2u
from numpy import alen, cumsum

def update_runtime_variables(model):
    if model.dynamic_h_init:
            model.h_init = min(model.c_gammah / model.gamma, model.h_init_max)
            if model.verbosity > 1:
                print('\nh_init: ', model.h_init,
                      'u_init', h2u(model.h_init, model.n,
                                    model.m, model.gamma), '\n')

def solve(model):

    def ip_direct_drainage(model, measurements_names):
        update_runtime_variables(model)
        (flag, t, z, gc1, rm1, u, wm, wm_in_tube) = solve_direct(model)

        contains_data = (alen(t) > 1)

        result = [flag, t]

        for name in measurements_names:
            # we discard values at t=0 (for given measurement)
            if name == 'MO':
                result.append(z[1:, model.mass_out_idx].transpose())
            elif name == 'MI':
                result.append(z[1:, model.mass_in_idx].transpose())
            elif name == 'GC':
                result.append(gc1[1:])
            elif name == 'RM':
                result.append(rm1[1:])

        return result

    t_meas = measurements_time(model)

    calc_p = model.get_parameters(('calc_gc', 'calc_rm', 'calc_wm')) # backup
    model.calc_gc = bool(model.gc1)
    model.calc_rm = bool(model.rm1)
    model.calc_wm = model.calc_gc or model.calc_rm

    (inv_params, cov) = \
      simulate_inverse(t_meas, ip_direct_drainage, model, model.inv_init_params,
                       wl_in_meas  = model.wl1,
                       wl_out_meas = model.wl_out,
                       gc_meas     = model.gc1,
                       rm_meas     = model.rm1,
                       wl_in_weights  = model.wl1_weights,
                       wl_out_weights = model.wl_out_weights,
                       gc_weights     = model.gc1_weights,
                       rm_weights     = model.rm1_weights,
                       optimfn=model.optimfn)

    model.set_parameters(calc_p) # ...and restore values

    return (inv_params, cov)

def run(model):
    (inv_params, cov) = solve(model)

    # DISPLAY RESULTS:
    if inv_params:
        model.set_parameters(inv_params)
        # run once again the direct problem with optimal parameters
        model.calc_wm = True
        model_verbosity = model.verbosity # backup verbosity
        model.verbosity = 0

        show_results(extract_data, model, inv_params=inv_params, cov=cov)

        model.verbosity = model_verbosity # restore verbosity

    return inv_params
