from sys import path as syspath
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from shared_functions import h2u

import modules.direct_draining_saturated.run as direct
from config import merge_cfgs

syspath.append('/'.join(['.', 'odes', 'build', 'lib.linux-x86_64-3.2']))

#n, gamma are dummy values, just to assure that direct checks will pass

PARAMETERS = {'inverse': ['inv_init_params']}
# we make sure that n is specified - otherwise direct.adjust fails
CFG_ADDITIONAL_PARAMETERS = {'soil': {'ks': 0.0, 'n': 1.0, 'gamma': 0.0}}


def base_cfg():
    return merge_cfgs(direct.base_cfg(), CFG_ADDITIONAL_PARAMETERS)

def adjust_cfg(flattened_cfg):
    direct.adjust_cfg(flattened_cfg)

    if not 'inv_init_params' in flattened_cfg:
        raise ValueError('CFG:check: Value not initialized: %s' 
                         % 'inv_init_params')
    if len(flattened_cfg['inv_init_params']) == 2 and flattened_cfg['ks'] == 0.0:
        raise ValueError('CFG:check: Value not initialized: %s' 
                         % 'Ks')
    #print(flattened_cfg)
    #raise ValueError('eee')

def solve(model):
    def ip_direct_drainage(model, optim_args):
        if len(optim_args) == 3:   # (Ks, n, gamma)
            (model.ks, model.n, model.gamma) = optim_args
        elif len(optim_args) == 2: # (n, gamma)
            (model.n, model.gamma) = optim_args
        model.m = 1-1/model.n

        (_flag, t, z) = direct.solve(model)

        u = h2u(z[:, model.first_idx:model.last_idx+1], 
                model.n, model.m, model.gamma)

        # characteristics: GC is measured only inside the tube, so no mass 
        #   on output is taken into accout
        mass_out = 0.0
        GC, _RM, _WM = \
          direct.characteristics(t, u,
                                 z[:, model.mass_in_idx], mass_out,
                                 z[:, model.s1_idx], z[:, model.s2_idx], model, 
                                 chtype='gc')

        print('\noptimization parameters: ', optim_args)
        wl = z[:, model.mass_out_idx]
        print('gc_mes, wl_mes: t_exp', model.gc1, model.wl_out1, model.duration)
        print('gc_com, wl_com: t_com', GC[1:], wl[1:], t)
        #input('pause...')

        return (t, GC[1:], wl.transpose()[1:]) # discard values at t=0

    def lsq_ip_direct_drainage(xdata, *optim_args):

        (t, GC, wl) = ip_direct_drainage(xdata, optim_args)
        #print('gc, wl: ', GC, wl)
        #print(np.shape(GC), np.shape(wl), type(GC), type(wl))
        result = np.concatenate((GC, wl))
        #print(result)
        return result


    #model.r0 = [model.r0_fall for wl0 in model.wl0]
    model.omega_fall = (np.sqrt(model.g/model.r0_fall)
                        * np.ones(np.shape(model.l0), float))

    # resolve the type of measured data
    if model.exp_type in ['ids', 'idsh']:
        #print(np.asarray(model.wl_out1), model.gc1)
        data_measured = np.concatenate((np.asarray(model.gc1, dtype=float),
                                        np.asarray(model.wl_out1, dtype=float)))
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

    inv_params_init = tuple(model.inv_init_params)

    # Solve inverse problem
    #    Ks_inv, cov_ks = curve_fit(lsq_direct_fn, xdata,
    #                           data_measured, p0 = Ks_init)
    inv_params, cov_ks = curve_fit(lsq_direct_fn, xdata,
                                   data_measured, p0 = inv_params_init)

    (t_inv, GC_inv, wl_inv) = direct_fn(xdata, inv_params)

    # Print results
    for i in np.arange(len(model.duration)):
        if model.wl_out1[i] == 0.0:
            wl_out_measured = 1.0e-10
        else:
            wl_out_measured = model.wl_out1[i]
        wl_out_computed = wl_inv[i]

        duration_measured =  model.duration[i]
        duration_computed = t_inv[i+1] - t_inv[i]
        gc_measured = model.gc1[i]
        gc_computed = GC_inv[i]

        print('Subexperiment ', i+1)
        print('    GC_measured: % 3.6f wl_out_measured: % .6f'
              ' t_end_expected: % 3.2f'
              % (gc_measured,  wl_out_measured, duration_measured))
        print('    GC_computed: % 3.6f wl_out_computed: % .6f'
              ' t_end_computed: % 3.2f'
              % (gc_computed, wl_out_computed, duration_computed))
        print('    Error (%%):   % 2.2f                        % 2.2f'
              '                  % 2.2f'
              % ((gc_computed - gc_measured) / gc_measured * 100,
                 (wl_out_computed - wl_out_measured) / wl_out_measured * 100,
                 (duration_computed - duration_measured) / duration_measured * 100))
    if inv_params_len == 3:
        (Ks_inv, n_inv, gamma_inv) = inv_params
        print('\nKs [cm/s] found: ', Ks_inv)
    elif inv_params_len == 2:
        (n_inv, gamma_inv) = inv_params
    print('n         found: ', n_inv)
    print('gamma     found: ', gamma_inv)

    return inv_params

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
