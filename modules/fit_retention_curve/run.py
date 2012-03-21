import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import modules.base.run as base
from config import merge_cfgs
from shared_functions import h2u

PARAMETERS = {'inverse': ['inv_init_params']}

CFG_ADDITIONAL_PARAMETERS = {}

def base_cfg():
    return CFG_ADDITIONAL_PARAMETERS

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

def experiments_files(first_experiment, last_experiment, tubes):
    files = []
    identifiers = []

    for exp_no in range(first_experiment, last_experiment+1):
            inifilename = 'experiment_' + str(exp_no) +'.ini'
            identifier  = 'experiment ' + str(exp_no)

            files.append(inifilename)
            identifiers.append(identifier)

    return (identifiers, files)


def solve(model):
    def lsq_fn(xdata, *optim_args):
        print(xdata, optim_args)

        optim_args_len = len(optim_args)

        if optim_args_len == 4:
            h = xdata
            print(h)
            print(optim_args)
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

        print('th', theta)

        return theta

    inv_params_init = model.inv_init_params
    inv_init_params_len = len(model.inv_init_params)

    if type(model.p) == list:
        h = np.asarray([-p / model.rho / model.g for p in model.p])
    else:
        h = model.p / model.rho /model.g

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

    print('theta_measured', model.theta)

    if inv_init_params_len == 4:
        (n, gamma, theta_s, theta_r) = inv_params
        print(' n:\t %s\n gamma:\t %s\n theta_s:\t %s\n theta_r:\t %s'
              % tuple(inv_params))
    elif inv_init_params_len == 3:
        if theta_s_p:
            (n, gamma, theta_r) = inv_params
            print(' n:\t %s\n gamma:\t %s\n theta_r:\t %s' % tuple(inv_params))
        else:
            (n, gamma, theta_s) = inv_params
            print(' n:\t %s\n gamma:\t %s\n theta_s:\t %s' % tuple(inv_params))
    else:
        print(' n:\t %s\n gamma:\t %s\n theta_s:\t \n theta_r:\t'
              % tuple(inv_params))
    print(' Cov: %s' % cov_inv)

    #theta_inv = lsq_fn(xdata, *inv_params)

    draw_graphs(1, model.p, model.theta, n, gamma, theta_s, theta_r,
                model.rho, model.g)

    return inv_params

def draw_graphs(fignum, p_measured,theta_measured,
                n, gamma, theta_s, theta_r, rho, g):
    import matplotlib.pyplot as plt

    plt.figure(fignum, figsize=(8, 4.5))

    p_calc = np.arange(0, 60000, 1)
    theta_calc = theta_r + (theta_s - theta_r) * h2u(-p_calc/rho/g,
                                                     n, 1.-1./n, gamma)

    plt.plot(theta_calc, p_calc, '-', theta_measured, p_measured, 'x',)
    plt.yscale('log')
    plt.ylabel('Pressure $p$ [Pa]')
    plt.xlabel('Water content ${\theta}$ ')

    plt.show()
