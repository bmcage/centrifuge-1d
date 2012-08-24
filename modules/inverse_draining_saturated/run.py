from modules.direct_draining_saturated.run import solve as solve_direct
from modules.shared.functions import measurements_time
from modules.shared.solver import simulate_inverse
from modules.shared.vangenuchten import h2u
from numpy import alen

def solve(model):

    no_measurements = []

    def ip_direct_drainage(model):
        if model.dynamic_h_init:
            model.h_init = min(model.c_gammah / model.gamma, model.h_init_max)
            if model.verbosity > 1:
                print('\nh_init: ', model.h_init,
                      'u_init', h2u(model.h_init, model.n,
                                    model.m, model.gamma), '\n')

        (flag, t, z, gc1, rm1) = solve_direct(model)

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

    inv_params = \
      simulate_inverse(t_meas, ip_direct_drainage, model, model.inv_init_params,
                       wl_in_meas  = model.get_iterable_value('wl1'),
                       wl_out_meas = model.get_iterable_value('wl_out'),
                       gc_meas     = model.get_iterable_value('gc1'),
                       rm_meas     = model.get_iterable_value('rm1'),
                       optimfn=model.optimfn)

    return inv_params

def run(model):
    return solve(model)
