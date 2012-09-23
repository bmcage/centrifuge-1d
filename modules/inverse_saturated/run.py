from modules.direct_saturated.run import solve as direct_solve
from modules.shared.solver import simulate_inverse
from numpy import alen

def solve(model):

    def ip_direct_saturated_heights(model, measurements_names):

        (flag, t, z) = direct_solve(model)

        result = [flag, t]

        for name in measurements_names:
            if name == 'MI':
                result.append(z[1:, model.mass_in_idx].transpose())
            elif name == 'MO':
                result.append(z[1:, model.mass_out_idx].transpose())

        return result

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
