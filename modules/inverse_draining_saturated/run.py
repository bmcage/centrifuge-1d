from modules.direct_draining_saturated.run import solve as solve_direct
from modules.shared.solver import simulate_inverse
from modules.shared.shared_functions import (scale_array,
                                             determine_scaling_factor)
from modules.shared.show import disp_inv_results
from numpy import concatenate, asarray, empty, log, exp, cumsum,

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


    def ip_direct_drainage(model, optim_args, retn_full_p = False):

        if determine_all:
            (log_ks, log_n, log_gamma) = optim_args
            ks = exp(log_ks) / model.ks_inv_scale
            model.ks = ks
        else:
            (log_n, log_gamma) = optim_args

        n     = 1+exp(log_n)
        gamma = -exp(log_gamma)

        model.n     = n
        model.m     = 1-1/n
        model.gamma = gamma

        if model.verbosity > 0:
            print()
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
            (flag, t, z, gc1, rm1) = solve_direct(model)

            if not flag:
                print('Solver could not continue. Exiting...')
                exit(1)

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

        if model.verbosity > 0:
            disp_inv_results(model, t, wl_out_inv=wl_out, gc1_inv=gc1,
                             n_inv=None, gamma_inv=None,ks_inv=None,
                             display_graphs=False)

        if retn_full_p:
            return (t, wl_out, gc1, rm1, z)
        else:
            return (t, wl_out, gc1, rm1)

    data_measured = concatenate((t_meas, wl_out_meas, gc_meas, rm_meas))
    xdata         = model

    if determine_all:
        (ks_init, n_init, gamma_init) = model.inv_init_params
        ks_init = ks_init * model.ks_inv_scale
        init_params = (log(ks_init), log(n_init - 1.0), log(-gamma_init))
    else:
        (n_init, gamma_init) = model.inv_init_params
        init_params  = (log(n_init - 1.0), log(-gamma_init))

    inv_params, cov_ks = simulate_inverse(ip_direct_drainage, xdata,
                                          data_measured, init_params,
                                          optimfn=model.optimfn)

    (t_inv, wl_out_inv, gc1_inv, rm1_inv, z) = \
      ip_direct_drainage(xdata, inv_params, retn_full_p = True)

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
