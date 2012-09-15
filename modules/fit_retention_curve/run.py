import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import modules.base.run as base
from modules.shared.vangenuchten import h2u
from modules.shared.show import make_dplot, add_dplotline, display_dplot

def solve(model):
    def lsq_fn(xdata, *optim_args):
        print(optim_args)

        optim_args_len = len(optim_args)

        if optim_args_len == 4:
            h = xdata
            (n, gamma, theta_s, theta_r) = optim_args
        elif optim_args_len == 3:
            (h, theta_s_p, theta_sr) = xdata
            if theta_s_p:
                theta_s = theta_sr
                (n, gamma, theta_r) = optim_args
            else:
                theta_r = theta_sr
                (n, gamma, theta_s) = optim_args
        else:
            (h, theta_s, theta_r) = xdata
            (n, gamma) = optim_args

        m = 1. - 1./n

        u = h2u(h, n, m, gamma)

        theta = theta_r + (theta_s - theta_r) * u

        return theta

    inv_params_init = model.inv_init_params
    inv_init_params_len = len(model.inv_init_params)

    # [p] = Pa = kg/m/s^2 = 10 * g/cm/s^2 -\
    # [h] = cm                            - \
    # => h [cm] =(10*p)/g/rho with [g]=cm/s^2, [rho]=g/cm^3

    if type(model.p) == list:
        h = np.asarray([-10.*p / model.rho / model.g for p in model.p])
    else:
        h = 10.*model.p / model.rho /model.g

    if inv_init_params_len == 4:
        xdata = h
    elif inv_init_params_len == 3:
        theta_s_p = hasattr(model, 'theta_s')
        if theta_s_p:
            theta_s  = model.theta_s
            theta_sr = theta_s
        else:
            theta_r  = model.theta_r
            theta_sr = theta_r

        xdata = (h, theta_s_p, theta_sr)
    else:
        theta_s = model.theta_s
        theta_r = model.theta_r

        xdata = (h, theta_s, theta_r)

    data_measured = model.theta

    inv_params, cov_inv = curve_fit(lsq_fn, xdata,
                                    data_measured, p0 = inv_params_init)

    theta_inv = lsq_fn(xdata, *inv_params)
    print('\ntheta_measured:', end='')
    for th_m in data_measured:
        print(' % 8.4f' % th_m, end='')
    print('\ntheta_computed:', end='')
    for th_c in theta_inv:
        print(' % 8.4f' % th_c, end='')
    print('\nError  (%):    ', end='')
    for (th_m, th_c) in zip(data_measured, theta_inv):
        print(' % 8.4f' % ((th_c - th_m) / th_m * 100), end='')

    print('\n\nOptimized parameters found:')
    if inv_init_params_len == 4:
        (n, gamma, theta_s, theta_r) = inv_params
        params = ['n', 'gamma', 'theta_s', 'theta_r']
    elif inv_init_params_len == 3:
        if theta_s_p:
            (n, gamma, theta_r) = inv_params
            params = ['n', 'gamma', 'theta_r']
        else:
            (n, gamma, theta_s) = inv_params
            params = ['n', 'gamma', 'theta_s']
    else:
        (n, gamma) = inv_params
        params = ['n', 'gamma']

    for (param, value) in zip(params, inv_params):
        if param == 'gamma':
            print(' %-7s: % 8.5g' % (param, value))
        else:
            print(' %-7s: % 8.5f' % (param, value))

    print('\n Cov:\n%s\n' % cov_inv)

    def compute_theta(p, n, m, gamma, theta_s, theta_r, rho,g):
        theta_calc = theta_r + (theta_s - theta_r) * h2u(-10.*p/rho/g,
                                                         n, 1.-1./n, gamma)
        return theta_calc

    if model.show_figures:
        p_calc = np.arange(0, 10000000, 100)

        dplot = make_dplot('RC', legend_loc=1, yscale='log')

        # computed data
        theta_calc = compute_theta(p_calc, n, 1-1/n, gamma,
                                   theta_s, theta_r, model.rho, model.g)
        add_dplotline(dplot, theta_calc, p_calc, line_opts='-',
                      label='computed')
        # measured data
        if model.p and model.theta:
            add_dplotline(dplot, model.theta, model.p, line_opts='x',
                          label='measured')

        # referencing data
        n_ref     = model.n_ref
        gamma_ref = model.gamma_ref
        if n_ref and gamma_ref:
            theta_ref = compute_theta(p_calc, n_ref, 1-1/n_ref, gamma_ref,
                                      theta_s, theta_r, model.rho, model.g)
            add_dplotline(dplot, theta_ref, p_calc, line_opts='-')

        display_dplot(dplot, show_figures=True, separate_figures=True,)

    return inv_params

def run(model):
    return solve(model)
