import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from modules.shared.vangenuchten import h2u

def solve(model):
    def lsq_fn(xdata, *optim_args):
        print(optim_args)

        optim_args_len = len(optim_args)

        if optim_args_len == 3:
            (h, theta_s) = xdata
            (n, gamma, theta_r) = optim_args
        else:
            (h, theta_s, theta_r) = xdata
            (n, gamma) = optim_args

        m = 1. - 1./n

        u = h2u(h, n, m, gamma)
        #print(u, theta_r, theta_s)
        theta = theta_r + (theta_s - theta_r) * u

        return theta

    inv_params_init = model.inv_init_params
    inv_init_params_len = len(model.inv_init_params)

    re = np.asarray(model.re)
    l1 = np.asarray(model.l1)

    wm0 = model.porosity[0] * model.l0[0]
    wm_in_tube = wm0 - np.asarray(model.wl_out1)
    V_total = l1
    theta = wm_in_tube / V_total

    h = (-np.power(np.asarray(model.omega), 2) / np.asarray(model.g) / 2
         * (np.power(re, 2) - np.power(re - l1/2, 2)))
    h[0] = -0.001
    print(h)

    theta_s  = np.asarray(model.porosity)

    if inv_init_params_len == 3:
        xdata = (h, theta_s)
    else:
        theta_r = model.theta_r

        xdata = (h, theta_s, theta_r)

    data_measured = theta
    print('theta:', theta, theta_s, theta_r)

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
    if inv_init_params_len == 3:
        (n, gamma, theta_r) = inv_params
        params = ['n', 'gamma', 'theta_r']
    else:
        (n, gamma) = inv_params
        params = ['n', 'gamma']

    for (param, value) in zip(params, inv_params):
        print(' %-7s: % 8.5f' % (param, value))

    print('\n Cov:\n%s\n' % cov_inv)

    if model.draw_graphs:
        draw_graphs(n, gamma, theta_s, theta_r,
                    model.rho, model.g,
                    p_measured=-h*model.rho*model.g,
                    theta_measured=theta)

    return inv_params

def draw_graphs(n, gamma, theta_s, theta_r, rho, g,
                p_measured = None, theta_measured = None, fignum = 1):
    import matplotlib.pyplot as plt

    plt.figure(fignum, figsize=(8, 4.5))
    theta_s = theta_s[0]

    p_calc = np.logspace(-10, 8, 100)
    theta_calc = theta_r + (theta_s - theta_r) * h2u(-10.*p_calc/rho/g,
                                                     n, 1.-1./n, gamma)

    plt.plot(theta_calc, p_calc, '-')
    if (not p_measured is None) and (not theta_measured is None):
        plt.plot(theta_measured, p_measured, 'x',)
    plt.yscale('log')
    plt.ylabel('Pressure $p$ [Pa]')
    plt.xlabel('Water content \t$\theta$ ')

    plt.show()
