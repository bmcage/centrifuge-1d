from modules.direct_draining_saturated.run import solve as solve_direct
from modules.shared.functions import measurements_time
from modules.shared.solver import simulate_inverse
from modules.shared.vangenuchten import h2u
from modules.shared.show import disp_inv_results
from numpy import alen

def update_runtime_variables(model):
    if model.dynamic_h_init:
            model.h_init = min(model.c_gammah / model.gamma, model.h_init_max)
            if model.verbosity > 1:
                print('\nh_init: ', model.h_init,
                      'u_init', h2u(model.h_init, model.n,
                                    model.m, model.gamma), '\n')

def solve(model):

    no_measurements = []

    def ip_direct_drainage(model):
        update_runtime_variables(model)
        (flag, t, z, gc1, rm1, u, wm, wm_in_tube) = solve_direct(model)

        contains_data = (alen(t) > 1)

        # we discard values at t=0 (for given measurement)
        if model.calc_wl_out and contains_data:
            wl_out = z[1:, model.mass_out_idx].transpose()
        else:
            wl_out = no_measurements

        if model.calc_wl_in and contains_data:
            wl_in = z[1:, model.mass_in_idx].transpose()
        else:
            wl_in = no_measurements

        if model.calc_gc and contains_data:
            gc1 = gc1[1:]
        else:
            gc1 = no_measurements

        if model.calc_rm and contains_data:
            rm1 = rm1[1:]
        else:
            rm1 = no_measurements

        return (flag, t, wl_in, wl_out, gc1, rm1)

    t_meas = measurements_time(model)

    (inv_params, cov) = \
      simulate_inverse(t_meas, ip_direct_drainage, model, model.inv_init_params,
                       wl_in_meas  = model.get_iterable_value('wl1'),
                       wl_out_meas = model.get_iterable_value('wl_out'),
                       gc_meas     = model.get_iterable_value('gc1'),
                       rm_meas     = model.get_iterable_value('rm1'),
                       wl_in_weights  = model.get_iterable_value('wl1_weights'),
                       wl_out_weights = model.get_iterable_value('wl_out_weights'),
                       gc_weights     = model.get_iterable_value('gc1_weights'),
                       rm_weights     = model.get_iterable_value('rm1_weights'),
                       optimfn=model.optimfn)

    return (inv_params, cov)

def run(model):
    (inv_params, cov) = solve(model)

    # DISPLAY RESULTS:
    if inv_params:
        # run once again the direct problem with optimal parameters
        model.calc_wm = True
        model_verbosity = model.verbosity # backup verbosity
        model.verbosity = 0
        (flag, t, z, gc1, rm1, u, wm, wm_in_tube) = solve_direct(model)
        model.verbosity = model_verbosity # restore verbosity

        disp_inv_results(model, t, inv_params=inv_params, cov=cov,
                         wl_in_inv=z[1:, model.mass_in_idx].transpose(),
                         wl_out_inv=z[1:, model.mass_out_idx].transpose(),
                         gc1_inv=gc1, rm1_inv=rm1, wm=wm, y=model.y,
                         h_inv=z[:, model.first_idx:model.last_idx+1], u_inv=u,
                         s1_inv=z[:, model.s1_idx], s2_inv=z[:, model.s2_idx],
                         disp_abserror=True, display_graphs=True)

    return inv_params
