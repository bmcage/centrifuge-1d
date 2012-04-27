from modules.direct_draining_saturated.run import solve as solve_direct
from modules.shared.solver import simulate_inverse
from modules.shared.shared_functions import scale_array
from numpy import concatenate, asarray, empty

def print_results(model, inv_params, t_inv, wl_out1_inv, gc1_inv):
    duration = model.get_iterable_value('duration')

    subexperiments = len(duration)

    if model.calc_wl_out:
        wl_out1  = np.asarray(model.get_iterable_value('wl_out1'))
        wl_out1[wl_out1 == 0.0] = 1.0e-10
        print('WL_out_measured: ', subexperiments *  '% 8.6f'
              % tuple(wl_out1))
        print('WL_out_computed: ', subexperiments *  '% 8.6f'
              % tuple(wl_out1_inv))
        print('Error (\%):  ', subexperiments * '    % 5.2f'
              % tuple((gc1_inv - gc1) / gc1 * 100.))
    if model.calc_gc:
        gc1  = np.asarray(model.get_iterable_value('gc1'), dtype=float)
        print('GC_measured: ', subexperiments *  '% 9.6f' % tuple(gc1))
        print('GC_computed: ', subexperiments *  '% 9.6f' % tuple(gc1_inv))
        print('Error (\%):  ', subexperiments * '    % 5.2f'
              % tuple((gc1_inv - gc1) / gc1 * 100.))

    if model.calc_rm:
        rm1  = np.asarray(model.get_iterable_value('rm1'), dtype=float)
        print('RM_measured: ', subexperiments *  '% 9.6f' % tuple(rm1))
        print('RM_computed: ', subexperiments *  '% 9.6f' % tuple(rm1_inv))
        print('Error (\%):  ', subexperiments * '    % 5.2f'
              % tuple((rm1_inv - rm1) / rm1 * 100.))

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

        (_flag, t, z, gc1, rm1) = solve_direct(model)

        # we discard values at t=0 (for give measurement)
        if model.calc_wl_out:
            wl_out1 = z[1:, model.mass_out_idx].transpose()
        else:
            wl_out1 = empty([0,], dtype=float)

        if model.calc_gc:
            gc1 = gc1[1:]
        if model.calc_rm:
            rm1 = rm1[1:]

        print('\noptimization parameters: ', optim_args)
        print('gc_mes, wl_mes: t_exp', model.get_iterable_value('gc1'),
              model.get_iterable_value('wl_out1'),
              model.get_iterable_value('duration'))
        print('gc_com, wl_com: t_com', gc1[1:], wl_out1[1:], t)

        return (t, wl_out1, gc1, rm1)

    def lsq_ip_direct_drainage(model, *optim_args):

        (t, wl_out1, gc1, rm1) = ip_direct_drainage(model, optim_args)

        if model.calc_wl_out:
            scale_array(wl_out1, wl_out1)
        if model.calc_gc:
            scale_array(gc1, gc1)
        if model.calc_rm:
            scale_array(rm1, rm1)

        result = concatenate((wl_out1, gc1, rm1))

        return result

    if model.exp_type in ['ids', 'idsh']:
        if model.calc_wl_out:
            wl_out_meas = asarray(model.get_iterable_value('wl_out1'),
                                  dtype=float)
            scale_array(wl_out_meas, wl_out_meas)
        else:
            wl_out_meas = empty([0,], dtype=float)

        if model.calc_gc:
            gc_meas  = asarray(model.get_iterable_value('gc1'), dtype=float)
            scale_array(gc_meas, gc_meas)
        else:
            gc_meas = empty([0,], dtype=float)

        if model.calc_rm:
            rm_meas  = asarray(model.get_iterable_value('rm1'), dtype=float)
            scale_array(rm_meas, rm_meas)
        else:
            rm_meas = empty([0,], dtype=float)


        data_measured = concatenate((wl_out_meas, gc_meas, rm_meas))

        lsq_direct_fn = lsq_ip_direct_drainage
        direct_fn     = ip_direct_drainage

    xdata         = model
    init_params   = model.inv_init_params

    inv_params, cov_ks = simulate_inverse(lsq_direct_fn, xdata, data_measured,
                                          init_params)

    (t_inv, wl_out1_inv, gc1_inv, rm1_inv) = direct_fn(xdata, inv_params)

    print_results(model, inv_params, t_inv, wl_out1_inv, gc1_inv)

    return inv_params
