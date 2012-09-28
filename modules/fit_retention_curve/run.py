import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import modules.base.run as base
from modules.shared.vangenuchten import h2u
from modules.shared.show import make_dplot, add_dplotline, display_dplots

def solve(model):

    def lsq_fn(xdata, *optim_args):
        print(optim_args)

        for (name, value) in zip(params_names, optim_args):
            xdata[name] = value

        n = xdata['n']
        m = 1. - 1./n

        u = h2u(xdata['h'], n, m, xdata['gamma'])

        (theta_s, theta_r) = (xdata['theta_s'], xdata['theta_r'])
        theta = theta_r + (theta_s - theta_r) * u

        return theta

    # [p] = Pa = kg/m/s^2 = 10 * g/cm/s^2 -\
    # [h] = cm                            - \
    # => h [cm] =(10*p)/g/rho with [g]=cm/s^2, [rho]=g/cm^3
    p = np.asarray(model.p, dtype=float)
    h = -10.*p / model.rho /model.g

    xdata = {'h': h, 'theta_s': 0.0, 'theta_r': 0.0, 'gamma': 0.0, 'n': 0.0}

    (params_names, inv_params_init) = ([], [])
    for name in ('theta_s', 'theta_r', 'n', 'gamma'):
        if name in model.inv_init_params:
            params_names.append(name)
            inv_params_init.append(model.inv_init_params[name])
        elif hasattr(model, name):
            xdata[name] = getattr(model, name)
        else:
            print('Unknown parameter: ', name, '\nExiting.')
            exit(1)

    data_measured = model.theta

    (inv_params, cov_inv) = \
      curve_fit(lsq_fn, xdata, data_measured, p0 = inv_params_init)

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

        dplot = make_dplot('RC', legend_loc=1, yscale='log', legend_title=None)

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
            add_dplotline(dplot, theta_ref, p_calc, line_opts='-',
                          label=model.label_ref)

        display_dplots(dplot, show_figures=True, separate_figures=True,)

    return inv_params

def run(model):
    return solve(model)
