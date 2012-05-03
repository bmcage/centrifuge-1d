from modules.direct_draining_saturated.run import solve as solve_direct
from modules.shared.solver import simulate_inverse
from modules.shared.shared_functions import (scale_array,
                                             determine_scaling_factor)
from numpy import concatenate, asarray, empty, zeros, ones, log, exp, cumsum

def print_results(model, wl_out1_inv, gc1_inv, t_inv, n_inv, gamma_inv,
                  ks_inv = None):
    duration = model.get_iterable_value('duration')

    subexperiments = len(duration)

    if model.calc_wl_out:
        wl_out1  = asarray(model.get_iterable_value('wl_out1'))
        wl_out1[wl_out1 == 0.0] = 1.0e-10
        print('WL_out_measured: ', subexperiments *  '% 8.6f'
              % tuple(wl_out1))
        print('WL_out_computed: ', subexperiments *  '% 8.6f'
              % tuple(wl_out1_inv))
        print('Error (\%):  ', subexperiments * '    % 5.2f'
              % tuple((wl_out1_inv - wl_out1) / wl_out1 * 100.))
    if model.calc_gc:
        gc1  = asarray(model.get_iterable_value('gc1'), dtype=float)
        print('GC_measured: ', subexperiments *  '% 9.6f' % tuple(gc1))
        print('GC_computed: ', subexperiments *  '% 9.6f' % tuple(gc1_inv))
        print('Error (\%):  ', subexperiments * '    % 5.2f'
              % tuple((gc1_inv - gc1) / gc1 * 100.))

    if model.calc_rm:
        rm1  = asarray(model.get_iterable_value('rm1'), dtype=float)
        print('RM_measured: ', subexperiments *  '% 9.6f' % tuple(rm1))
        print('RM_computed: ', subexperiments *  '% 9.6f' % tuple(rm1_inv))
        print('Error (\%):  ', subexperiments * '    % 5.2f'
              % tuple((rm1_inv - rm1) / rm1 * 100.))

    if not ks_inv is None:
        print('\nKs [cm/s] found: ', ks_inv)
    print('n         found: ', n_inv)
    print('gamma     found: ', gamma_inv)

def solve(model):
    # use all 3 parameters (Ks, n, gamma) or only two (n, gamma) ?
    determine_all = (len(model.inv_init_params) == 3)

    # prepare the measured data
    if model.calc_wl_out:
        wl_out_meas = asarray(model.get_iterable_value('wl_out1'), dtype=float)
        c_coef_wl_out = determine_scaling_factor(wl_out_meas)
        scale_array(wl_out_meas, c_coef_wl_out, wl_out_meas)
    else:
        wl_out_meas = empty([0,], dtype=float)

    if model.calc_gc:
        gc_meas  = asarray(model.get_iterable_value('gc1'), dtype=float)
        c_coef_gc = determine_scaling_factor(gc_meas)
        scale_array(gc_meas, c_coef_gc, gc_meas)
    else:
        gc_meas = empty([0,], dtype=float)

    if model.calc_rm:
        rm_meas  = asarray(model.get_iterable_value('rm1'), dtype=float)
        c_coef_rm = determine_scaling_factor(rm_meas)
        scale_array(rm_meas, c_coef_rm, rm_meas)
    else:
        rm_meas = empty([0,], dtype=float)

    t_meas = cumsum([0] + model.get_iterable_value('duration'), dtype=float)



    def ip_direct_drainage(model, optim_args):
        if determine_all:
            (model.ks, log_n, log_gamma) = optim_args
        else:
            (log_n, log_gamma) = optim_args

        model.n     = 1+exp(log_n)
        model.m     = 1-1/model.n
        model.gamma = -exp(log_gamma)

        if determine_all:
            print('Ks:    ', model.ks)
        print('n:     ', model.n)
        print('gamma:', model.gamma)
        # input('ENTER...')

        if (model.n > 20.) or (model.gamma > -1e-8) or (model.gamma < -4):
            # untolerable range, the solver will probably crash so we
            # return values based on how far from expected range we are
            if model.gamma > -1e-8:
                gamma_factor = -1./(1e8 * model.gamma)
            else:
                gamma_factor = exp(model.gamma)

            n_factor     = exp(model.n)

            if model.calc_wl_out:
                wl_out1 = (model.get_iterable_value('wl_out1')
                           + n_factor + gamma_factor)
            else:
                wl_out1 = empty([0,], dtype=float)

            if model.calc_gc:
                gc1 = (model.get_iterable_value('gc1')
                           + n_factor + gamma_factor)
            else:
                gc1 = empty([0,], dtype=float)

            if model.calc_rm:
                rm1 = (model.get_iterable_value('rm1')
                           + n_factor + gamma_factor)
            else:
                rm1 = empty([0,], dtype=float)

            t = t_meas
        else:
            (_flag, t, z, gc1, rm1) = solve_direct(model)

            # we discard values at t=0 (for give measurement)
            if model.calc_wl_out:
                wl_out1 = z[1:, model.mass_out_idx].transpose()
                scale_array(wl_out1, c_coef_wl_out, wl_out1)
            else:
                wl_out1 = empty([0,], dtype=float)

            if model.calc_gc:
                gc1 = gc1[1:]
                scale_array(gc1, c_coef_gc, gc1)
            if model.calc_rm:
                rm1 = rm1[1:]
                scale_array(rm1, c_coef_rm, rm1)

        print('gc_mes, wl_mes, t_exp', gc_meas, wl_out_meas, t_meas)
        print('gc_com, wl_com, t_com', gc1, wl_out1, t)

        return (t, wl_out1, gc1, rm1)

    data_measured = concatenate((t_meas, wl_out_meas, gc_meas, rm_meas))
    xdata         = model

    if determine_all:
        (ks_init, n_init, gamma_init) = model.inv_init_params
        init_params = (ks_init, log(n_init - 1.0), log(-gamma_init))
    else:
        (n_init, gamma_init) = optim_args
        init_params  = (log(n_init - 1.0), log(-gamma_init))

    inv_params, cov_ks = simulate_inverse(ip_direct_drainage, xdata,
                                          data_measured, init_params)

    (t_inv, wl_out1_inv, gc1_inv, rm1_inv) = ip_direct_drainage(xdata,
                                                                inv_params)

    if determine_all:
        (ks_inv, log_n, log_gamma) = inv_params
    else:
        (log_n, log_gamma) = inv_params
        ks_inv = None
    n_inv     = 1+exp(log_n)
    gamma_inv = -exp(log_gamma)

    if determine_all:
        resulting_params = (ks_inv, n_inv, gamma_inv)
    else:
        resulting_params = (n_inv, gamma_inv)

    if model.calc_wl_out:
        scale_array(wl_out1_inv, 1/c_coef_wl_out, wl_out1_inv)
    if model.calc_gc:
        scale_array(gc1_inv, 1/c_coef_gc, gc1_inv)
    if model.calc_rm:
        scale_array(rm1_inv, 1/c_coef_gc, rm1_inv)

    print_results(model, wl_out1_inv, gc1_inv, t_inv, n_inv, gamma_inv, ks_inv)

    return resulting_params
