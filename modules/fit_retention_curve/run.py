import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import modules.base.run as base
from shared_functions import h2u

PARAMETERS = {'inverse': ['inv_init_params']}

def adjust_cfg(flattened_cfg):

    required_parameters = ['inv_init_params', 'p', 'theta', 'draw_graphs',
                           'rho', 'g']

    for param in required_parameters:
        if not param  in flattened_cfg:
            print('CFG:check: Missing paramter in configuration file(s): %s'
                  % param)
            exit(1)

    theta_s_p = 'theta_s' in flattened_cfg
    theta_r_p = 'theta_r' in flattened_cfg
    inv_init_params_len = len(flattened_cfg['inv_init_params'])

    if ((theta_s_p and theta_r_p and inv_init_params_len == 2)
        or (theta_s_p and (not theta_r_p) and inv_init_params_len == 3)
        or ((not theta_s_p) and theta_r_p and inv_init_params_len == 3)
        or ((not theta_s_p) and (not theta_r_p) and inv_init_params_len == 4)):

        # check correctness of input:
        # if 2 initial guesses given, optimize only n, gamma
        # if 3 given, we optimize also either theta_s or theta_r
        # if 4 we optimize also theta_s and theta_r
        # theta_s and theta_r;
        pass
    else:
        th_s_str = ''
        th_r_str = ''
        if theta_s_p: th_s_str = ' = ' + str(flattened_cfg['theta_s'])
        if theta_r_p: th_s_str = ' = ' + str(flattened_cfg['theta_r'])

        print('Inconsistent initial guesses inv_init_params = %s'
              % flattened_cfg['inv_init_params'])
        print("with 'theta_s'%s and 'theta_r'%s" % (th_s_str, th_r_str))
        exit(1)

    flattened_cfg['theta_s_p'] = theta_s_p

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
        print(' %-7s: % 8.5f' % (param, value))

    print('\n Cov:\n%s\n' % cov_inv)

    if model.draw_graphs:
        draw_graphs(model.p, model.theta, n, gamma, theta_s, theta_r,
                    model.rho, model.g)

    return inv_params

def draw_graphs(fignum, p_measured,theta_measured,
                n, gamma, theta_s, theta_r, rho, g, fignum = 1):
    import matplotlib.pyplot as plt

    plt.figure(fignum, figsize=(8, 4.5))

    p_calc = np.arange(0, 10000000, 100)
    theta_calc = theta_r + (theta_s - theta_r) * h2u(-10.*p_calc/rho/g,
                                                     n, 1.-1./n, gamma)

    plt.plot(theta_calc, p_calc, '-', theta_measured, p_measured, 'x',)
    plt.yscale('log')
    plt.ylabel('Pressure $p$ [Pa]')
    plt.xlabel('Water content ${\theta}$ ')

    plt.show()
