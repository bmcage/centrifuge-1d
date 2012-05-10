from modules.direct_draining_saturated.run import solve as solve_direct
from modules.shared.solver import simulate_inverse
from modules.shared.shared_functions import (scale_array,
                                             determine_scaling_factor)
from numpy import concatenate, asarray, empty, zeros, ones, log, exp, cumsum

def print_results(model, wl_out_inv, gc1_inv, t_inv, n_inv, gamma_inv,
                  ks_inv = None, display_graphs=True, fignum = 1):
    duration = model.get_iterable_value('duration')

    subexperiments = len(duration)

    if model.calc_wl_out:
        wl_out  = asarray(model.get_iterable_value('wl_out'))
        wl_out  = wl_out.cumsum()
        wl_out[wl_out == 0.0] = 1.0e-10
        print('WL_out_measured: ', subexperiments *  '% 8.6f'
              % tuple(wl_out))
        print('WL_out_computed: ', subexperiments *  '% 8.6f'
              % tuple(wl_out_inv))
        print('Error (\%):  ', subexperiments * '    % 5.2f'
              % tuple((wl_out_inv - wl_out) / wl_out * 100.))
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

    if display_graphs:
        import matplotlib.pyplot as plt

        def add_legend(lines):
            legend_data = ['measured', 'computed']

            plt.figlegend(h_lines, legend_data, 1, borderaxespad=0.0,
                          prop={'family': 'monospace'})

        t_duration = model.get_iterable_value('duration')
        if t_duration is None:
            t_duration = model.duration
        t_fh_duration = model.get_iterable_value('fh_duration')
        if t_fh_duration is None:
            t_fh_duration = model.fh_duration

        t_minutes = cumsum((asarray(t_duration, dtype=float)
                            + asarray(t_fh_duration, dtype=float))
                            / 60)
        t_inv_minutes = t_inv[1:] / 60

        figures_count = (int(model.calc_wl_out) + int(model.calc_gc)
                         + int(model.calc_gc))

        xlabel = ('Time [min]')

        if figures_count == 1:
            plt.figure(fignum)
            rows = 1
            cols = 1

        elif figures_count == 2:
            plt.figure(fignum, figsize=(16, 4.5))
            rows = 1
            cols = 2
        else:
            plt.figure(fignum, figsize=(16, 8.5))
            rows = 2
            cols = 2

        plt.subplots_adjust(wspace=0.15, left=0.06, right=0.85)

        current_column = 1

        if model.calc_wl_out:
            plt.subplot(rows, cols, current_column)

            h_lines = plt.plot(t_minutes, wl_out, '.',
                               t_inv_minutes, wl_out_inv, 'x')
            plt.xlabel(xlabel)
            plt.ylabel('Expelled water [cm]')

            add_legend(h_lines)

            current_column = current_column + 1

        if model.calc_gc:
            plt.subplot(rows, cols, current_column)

            h_lines = plt.plot(t_minutes, gc1, '.', t_inv_minutes, gc1_inv, 'x')
            plt.xlabel(xlabel)
            plt.ylabel('Gravitational centre [cm]')

            add_legend(h_lines)

            current_column = current_column + 1

        if model.calc_rm:
            if rows == 2:
                current_column = 1

            plt.subplot(rows, cols, current_column)

            h_lines = plt.plot(t_minutes, rm1, '.', t_inv_minutes, rm1_inv, 'x')
            plt.xlabel(xlabel)
            plt.ylabel('Rotational momentum [cm]')

            add_legend(h_lines)

        plt.show(block=False)
        input('Enter...')

def solve(model):
    # use all 3 parameters (Ks, n, gamma) or only two (n, gamma) ?
    determine_all = (len(model.inv_init_params) == 3)

    # prepare the measured data
    if model.calc_wl_out:
        wl_out_sep  = asarray(model.get_iterable_value('wl_out'), dtype=float)
        wl_out_meas = wl_out_sep.cumsum()
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
            (log_ks, log_n, log_gamma) = optim_args
            ks = exp(log_ks)
            model.ks = ks
        else:
            (log_n, log_gamma) = optim_args

        n     = 1+exp(log_n)
        gamma = -exp(log_gamma)

        model.n     = n
        model.m     = 1-1/n
        model.gamma = gamma

        if determine_all:
            print('Ks:    ', model.ks)
        print('n:     ', model.n)
        print('gamma:', model.gamma)
        # input('ENTER...')

        ubounds = model.params_ubounds
        lbounds = model.params_lbounds

        out_of_range = False

        if n > ubounds['n']:
            out_of_range = True
            n_factor     = 10*exp(n-ubounds['n'])
        elif n < lbounds['n']:
            out_of_range = True
            n_factor = 1/(n-1)
        else:
            n_factor = 0.

        if gamma > ubounds['gamma']:
            out_of_range = True
            gamma_factor = 10 * (ubounds['gamma'] / gamma)
        elif gamma < lbounds['gamma']:
            out_of_range = True
            gamma_factor = 10*exp(lbounds['gamma'] - gamma)
        else:
            gamma_factor = 0.

        if determine_all:
            if ks > ubounds['ks']:
                out_of_range = True
                ks_factor = 10 * exp(ks - ubounds['ks'])
            elif ks < lbounds['ks']:
                out_of_range = True
                if ks == 0.0:
                    ks = 1e-20
                ks_factor = 10 * (lbounds['ks'] / ks)
            else:
                ks_factor = 0.
        else:
            ks_factor = 0.

        if out_of_range:
            # untolerable range, the solver will probably crash so we
            # return values based on how far from expected range we are

            penalization = n_factor + gamma_factor + ks_factor

            if model.calc_wl_out:
                wl_out = wl_out_meas + penalization
            else:
                wl_out = empty([0,], dtype=float)

            if model.calc_gc:
                gc1 = gc_meas + penalization
            else:
                gc1 = empty([0,], dtype=float)

            if model.calc_rm:
                rm1 = rm_meas + penalization
            else:
                rm1 = empty([0,], dtype=float)

            t = t_meas
        else:
            (_flag, t, z, gc1, rm1) = solve_direct(model)

            # we discard values at t=0 (for give measurement)
            if model.calc_wl_out:
                wl_out = z[1:, model.mass_out_idx].transpose()
                scale_array(wl_out, c_coef_wl_out, wl_out)
            else:
                wl_out = empty([0,], dtype=float)

            if model.calc_gc:
                gc1 = gc1[1:]
                scale_array(gc1, c_coef_gc, gc1)
            if model.calc_rm:
                rm1 = rm1[1:]
                scale_array(rm1, c_coef_rm, rm1)

        print('gc_mes, wl_mes, t_exp', gc_meas, wl_out_meas, t_meas)
        print('gc_com, wl_com, t_com', gc1, wl_out, t)

        return (t, wl_out, gc1, rm1)

    data_measured = concatenate((t_meas, wl_out_meas, gc_meas, rm_meas))
    xdata         = model

    if determine_all:
        (ks_init, n_init, gamma_init) = model.inv_init_params
        init_params = (log(ks_init), log(n_init - 1.0), log(-gamma_init))
    else:
        (n_init, gamma_init) = optim_args
        init_params  = (log(n_init - 1.0), log(-gamma_init))

    inv_params, cov_ks = simulate_inverse(ip_direct_drainage, xdata,
                                          data_measured, init_params)

    (t_inv, wl_out_inv, gc1_inv, rm1_inv) = ip_direct_drainage(xdata,
                                                                inv_params)

    if determine_all:
        (log_ks, log_n, log_gamma) = inv_params
        ks_inv = exp(log_ks)
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
        scale_array(wl_out_inv, 1/c_coef_wl_out, wl_out_inv)
    if model.calc_gc:
        scale_array(gc1_inv, 1/c_coef_gc, gc1_inv)
    if model.calc_rm:
        scale_array(rm1_inv, 1/c_coef_gc, rm1_inv)

    print_results(model, wl_out_inv, gc1_inv, t_inv, n_inv, gamma_inv, ks_inv)

    return resulting_params
