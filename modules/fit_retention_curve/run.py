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
    parameters = ['inv_init_params', 'h', 'u']

    for param in parameters:
        if not param  in flattened_cfg:
            print('CFG:check: Missing paramter in configuration file(s): %s' 
                  % param)
            exit(1)

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
    def lsq_fn(h, *optim_args):

        (n, gamma) = optim_args
        m = 1. - 1./n

        u = h2u(h, n, m, gamma)

        return u

    xdata = model.h
    inv_params_init = tuple(model.inv_init_params)
    data_measured   = model.u

    inv_params, cov_inv = curve_fit(lsq_fn, xdata,
                                   data_measured, p0 = inv_params_init)

    print(' n     found: %s\n gamma found: %s' % tuple(inv_params))
    print(' Cov: %s' % cov_inv)

    u_inv = lsq_fn(xdata, *inv_params)

    draw_graphs(1, model.h, model.u, inv_params)

    return inv_params

def draw_graphs(fignum, h_measured, u_measured, inv_found):
    import matplotlib.pyplot as plt

    plt.figure(fignum, figsize=(8, 4.5))

    [n_found, gamma_found] = inv_found
    h_calc = -np.arange(0, 1600, 1)
    u_calc = h2u(h_calc, n_found, 1.-1./n_found, gamma_found)

    plt.plot(u_calc, -h_calc, '-', u_measured, -h_measured, 'x',)
    plt.yscale('log')
    plt.ylabel('Hydraulic pressure ''h'' [$Pa$]')
    plt.xlabel('Relative saturation ''u'' ')

    plt.show()
