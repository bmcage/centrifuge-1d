import numpy as np
from modules.direct_saturated.run import solve as direct_solve
from modules.shared.solver import simulate_inverse

def print_results(model, Ks_inv, t_inv, wl1_inv):
    wl0      = model.get_iterable_value('wl0')
    wl1      = model.get_iterable_value('wl1')
    if model.duration == 0.0:
        duration = model.get_iterable_value('fh_duration')
    else:
        duration = model.get_iterable_value('duration')

    for i in range(model.iterations):
        print('Subexperiment %i:' % (i+1))
        print('    wl0         : % 9.6f' % wl0[i])
        print('    wl1_measured: % 9.6f    t_end_expected: % 9.2f' %
              (wl1[i], duration[i]))
        print('    wl1_computed: % 9.6f    t_end_computed: % 9.2f' %
              (wl1_inv[i], t_inv[i]))
        print('    Error (%%)   :  % 5.2f                            % 5.2f' %
              ((wl1_inv[i] - wl1[i]) / wl1[i] * 100,
               (t_inv[i] - duration[i]) / duration[i] * 100))

    print('\nKs found: ', Ks_inv)

def solve(model):
    def ip_direct_saturated_heights(model, Ks):
        model.ks = Ks

        (_flag, t, z) = direct_solve(model)

        return t, z[1:, model.mass_in_idx]

    def lsq_ip_direct_saturated_heights(xdata, Ks):
        return ip_direct_saturated_heights(xdata, Ks)[1] # return only wl

    # resolve the type of measured data
    if model.exp_type in ['ish', 'ish-sc', 'ish-f']:
        data_measured = model.get_iterable_value('wl1')
        lsq_direct_fn = lsq_ip_direct_saturated_heights
        direct_fn     = ip_direct_saturated_heights
    else:
        raise ValueError('Unrecognized experiment type exp_type: %s'
                         % model.exp_type)
    xdata         = model

    Ks_inv, cov_ks = simulate_inverse(lsq_direct_fn, xdata, data_measured,
                                      model.inv_init_params)

    t_inv, wl1_inv = direct_fn(xdata, Ks_inv)

    print_results(model, Ks_inv, t_inv[1:], wl1_inv)

    return Ks_inv
