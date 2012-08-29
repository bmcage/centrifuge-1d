from modules.direct_saturated.run import solve as direct_solve
from modules.shared.solver import simulate_inverse
from numpy import alen

def solve(model):

    no_measurements = []

    def ip_direct_saturated_heights(model):

        (flag, t, z) = direct_solve(model)

        contains_data = (alen(t) > 1)

        if model.calc_wl_out and contains_data:
            wl_out = z[1:, model.mass_out_idx].transpose()
        else:
            wl_out = no_measurements

        if model.calc_wl_in and contains_data:
            wl_in = z[1:, model.mass_in_idx].transpose()
        else:
            wl_in = no_measurements

        return (flag, t, wl_in, wl_out)

    t_meas = measurements_time(model)

    wl1_meas = model.get_iterable_value('wl1')

    inv_params = \
      simulate_inverse(t_meas, ip_direct_saturated_heights,
                       model, model.inv_init_params,
                       wl_in_meas  = model.get_iterable_value('wl1'),
                       wl_out_meas = model.get_iterable_value('wl_out'),
                       optimfn=model.optimfn)

    return inv_params

def run(model):
    return solve(model)
