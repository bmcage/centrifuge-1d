from modules.direct_draining_saturated.run import solve as solve_direct
from modules.shared.solver import simulate_inverse
from modules.shared.shared_functions import scale_array
from numpy import concatenate, asarray

def print_results(model, inv_params, t_inv, wl1_inv, gc1_inv):
    duration = model.get_iterable_value('duration')
    gc1  = model.get_iterable_value('gc1')
    wl_out1  = model.get_iterable_value('wl_out1')

    for i in range(model.iterations):
        if wl_out1[i] == 0.0:
            wl_out1[i] = 1.0e-10

        print('Subexperiment %i:' % (i+1))
        print('    GC_measured: % 9.6f wl_out_measured: % 8.6f'
              ' t_end_expected: % 6.2f'
              % (gc1[i],  wl_out1[i], duration[i]))
        print('    GC_computed: % 9.6f wl_out_computed: % 8.6f'
              ' t_end_computed: % 6.2f'
              % (gc1_inv[i], wl_out1_inv, t_inv[i]))
        print('    Error (%%):   % 5.2f                        % 5.2f'
              '                  % 5.2f'
              % ((gc1_inv[i] - gc1[i]) / gc1[i] * 100,
                 (wl_out1[i] - wl_out1_inv) / wl_out1[i] * 100,
                 (t_inv[i] - t[i]) / t[i] * 100))

    if len(inv_params) == 3:
        (Ks_inv, n_inv, gamma_inv) = inv_params
        print('\nKs [cm/s] found: ', Ks_inv)
    else:
        (n_inv, gamma_inv) = inv_params
    print('n         found: ', n_inv)
    print('gamma     found: ', gamma_inv)

def solve(model):
    # use all 3 parameters (Ks, n, gamma) or only two (n, gamma) ?
    determine_all = (len(model.inv_init_params) == 3)

    def ip_direct_drainage(model, optim_args):
        if determine_all:
            (model.ks, model.n, model.gamma) = optim_args
        else:
            (model.n, model.gamma) = optim_args
        model.m = 1-1/model.n

        (_flag, t, z, gc1) = solve_direct(model)

        wl_out1 = z[:, model.mass_out_idx]

        print('\noptimization parameters: ', optim_args)
        print('gc_mes, wl_mes: t_exp', model.get_iterable_value('gc1'),
              model.get_iterable_value('wl_out1'),
              model.get_iterable_value('duration'))
        print('gc_com, wl_com: t_com', gc1[1:], wl_out1[1:], t)

        return (t, gc1[1:], wl_out1.transpose()[1:]) # discard values at t=0

    def lsq_ip_direct_drainage(xdata, *optim_args):

        (t, gc1, wl_out1) = ip_direct_drainage(xdata, optim_args)

        scale_array(gc1, gc1)
        scale_array(wl_out1, wl_out1)

        result = concatenate((gc1, wl_out1))

        return result

    if model.exp_type in ['ids', 'idsh']:
        gc_meas      = asarray(model.get_iterable_value('gc1'), dtype=float)
        wl_out_meas = asarray(model.get_iterable_value('wl_out1'), dtype=float)

        scale_array(gc_meas, gc_meas)
        scale_array(wl_out_meas, wl_out_meas)

        data_measured = concatenate((gc_meas, wl_out_meas))

        lsq_direct_fn = lsq_ip_direct_drainage
        direct_fn     = ip_direct_drainage

    xdata         = model
    init_params   = model.inv_init_params

    inv_params, cov_ks = simulate_inverse(lsq_direct_fn, xdata, data_measured,
                                          init_params)

    (t_inv, gc1_inv, wl_out1_inv) = direct_fn(xdata, inv_params)

    print_results(model, inv_params, t_inv, wl_out1_inv, gc1_inv)

    return inv_params
