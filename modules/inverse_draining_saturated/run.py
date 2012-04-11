import modules.direct_draining_saturated.run.solve as solve_direct
from numpy import concatenate

def print_results(model, inv_params, t_inv, wl1_inv, gc1_inv):
    t    = model.duration
    gc1  = model.get_iterable_value('gc1')
    wl1  = model.get_iterable_value('wl1')

    for i in range(model.iterations):
        if wl1[i] == 0.0:
            wl1[i] = 1.0e-10

        print('Subexperiment %i:' % i+1)
        print('    GC_measured: % 9.6f wl_out_measured: % 8.6f'
              ' t_end_expected: % 6.2f'
              % (gc1[i],  wl1[i], duration[i]))
        print('    GC_computed: % 9.6f wl_out_computed: % 8.6f'
              ' t_end_computed: % 6.2f'
              % (gc1_inv[i], wl1_inv, t_inv[i]))
        print('    Error (%%):   % 5.2f                        % 5.2f'
              '                  % 5.2f'
              % ((gc1_inv[i] - gc1[i]) / gc1[i] * 100,
                 (wl1[i] - wl1_inv) / wl1[i] * 100,
                 (t_inv[i] - t[i]) / t[i] * 100))

    if len(inv_params) == 3:
        (Ks_inv, n_inv, gamma_inv) = inv_params
        print('\nKs [cm/s] found: ', Ks_inv)
    else:
        (n_inv, gamma_inv) = inv_params
    print('n         found: ', n_inv)
    print('gamma     found: ', gamma_inv)

def solve(model):
    determine_all = (len(optim_args) == 3)   # (Ks, n, gamma) vs. (n, gamma)

    def ip_direct_drainage(model, optim_args):
        if determine_all:
            (model.ks, model.n, model.gamma) = optim_args
        else:
            (model.n, model.gamma) = optim_args
        model.m = 1-1/model.n

        (_flag, t, z, gc1) = direct.solve(model)

        wl1 = z[:, model.mass_out_idx]

        print('\noptimization parameters: ', optim_args)
        print('gc_mes, wl_mes: t_exp', model.gc1, model.wl_out1, model.duration)
        print('gc_com, wl_com: t_com', gc1[1:], wl[1:], t)

        return (t, gc1[1:], wl1.transpose()[1:]) # discard values at t=0

    def lsq_ip_direct_drainage(xdata, *optim_args):

        (t, gc1, wl1) = ip_direct_drainage(xdata, optim_args)

        result = concatenate((gc1, wl1))

        return result

    if model.exp_type in ['ids', 'idsh']:
        data_measured = \
          np.concatenate((np.asarray(model.get_iterable_value('gc1'),
                                     dtype=float),
                          np.asarray(model.get_iterable_value('wl1'),
                                     dtype=float)))

        lsq_direct_fn = lsq_ip_direct_drainage
        direct_fn     = ip_direct_drainage

    xdata         = model

    inv_params, cov_ks = simulate_inverse(lsq_direct_fn)

    (t_inv, GC_inv, wl1_inv) = direct_fn(xdata, inv_params)

    print_results(model, inv_params, t_inv, wl1_inv, gc1_inv)

    return inv_params
