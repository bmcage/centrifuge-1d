from sys import path as syspath
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

import modules.direct_draining_saturated.run as direct
from config import merge_cfgs

syspath.append('/'.join(['.', 'odes', 'build', 'lib.linux-x86_64-3.2']))

#n, gamma are dummy values, just to assure that direct checks will pass
CFG_ADDITIONAL_PARAMETERS = {}


def base_cfg():
    return merge_cfgs(direct.base_cfg(), CFG_ADDITIONAL_PARAMETERS)

def adjust_cfg(flattened_cfg):
    if not 'n' in flattened_cfg:
        # we make sure that n is specified - otherwise direct.adjust fails
        flattened_cfg['n'] = -1.0
    direct.adjust_cfg(flattened_cfg)

    print(flattened_cfg)
    raise ValueError('eee')

def solve(model):
    def ip_direct_drainage(model, optim_args):
        if len(optim_args) == 3:   # (Ks, n, gamma)
            (model.ks, model.n, model.gamma) = optim_args
        elif len(optim_args) == 2: # (n, gamma)
            (model.n, model.gamma) = optim_args

        (_flag, t, z) = direct.solve(model)

        u = h2u(z[:, first_idx:last_idx+1], model.n, model.m, model.gamma)

        GC, _RM, _WM = \
          direct.characteristics(t, u,
                                 z[:, mass_in_idx], z[:, mass_out_idx],
                                 z[:, s1_idx], z[:, s2_idx], model, 
                                 chtype='gc')
        

        return (t, GC)

    def lsq_ip_direct_drainage(xdata, optim_args):

        (t, GC) = ip_direct_drainage(xdata, optim_args)
        return GC


    #model.r0 = [model.r0_fall for wl0 in model.wl0]
    model.omega_fall = (np.sqrt(model.g/model.r0_fall)
                        * np.ones(np.shape(model.l0), float))

    # resolve the type of measured data
    if exp_type in ['ids', 'idsh']:
        data_measured = model.gc
        lsq_direct_fn = lsq_ip_direct_drainage
        direct_fn     = ip_direct_drainage
    else:
        raise ValueError('Unrecognized experiment type exp_type: %s'
                         % model.exp_type)
    xdata         = model

    inv_params_len = len(model.inv_init_params)
    if inv_params_len == 3:
        [model.ks, model.n, model.gamma]
    elif inv_params_len == 2:
        [model.n, model.gamma]
    else:
        raise ValueError('InitParams: should be of length 3 or 2.')
    model.m = 1.0 - 1.0/model.n

    inv_params_init = tuple(

    # Solve inverse problem
    #    Ks_inv, cov_ks = curve_fit(lsq_direct_fn, xdata,
    #                           data_measured, p0 = Ks_init)
    params_inv, cov_ks = curve_fit(lsq_direct_fn, xdata,
                                   data_measured, p0 = [Ks_init])

    t_inv, GC_inv = direct_fn(xdata, Ks_inv)

    # Print results
    for i in np.arange(len(data_measured)):
        print('Subexperiment ', i+1)
        print('    GC_measured: % .6f    t_end_expected: % .2f' %
              (model.gc[i],  model.duration[i]))
        print('    GC_computed: % .6f    t_end_computed: % .2f' %
              (GC_inv[i], t_inv[i]))
        print('    Error (%%):   % .2f                        % .2f' %
              ((GC_inv[i] - model.gc[i]) / model.gc[i] * 100,
               (t_inv[i] - model.duration[i]) / model.duration[i] * 100))
    if len(init_params) == 3:
        (Ks_inv, n_inv, gamma_inv) = params_inv
        print('\nKs    found: ', Ks_inv)
    elif len(init_params) == 2:
        (n_inv, gamma_inv) = params_inv
    print('\nn     found: ', n_inv)
    print('\ngamma found: ', gamma_inv)

    

    return params_inv

def verify_inverse_data(model):
    if not model.inverse_data_filename:
        raise ValueError('Data file for inverse problem not specified !')

    if isinstance(model.inverse_data_filename, str):
        data_filenames = [model.inverse_data_filename]
    elif isinstance(model.inverse_data_filename, (list, tuple)):
        data_filenames = model.inverse_data_filename
    else:
        raise ValueError('Wrong inverse_data_filename: %d'
                         % model.data_type)

if __name__ == "__main__":
    pass
