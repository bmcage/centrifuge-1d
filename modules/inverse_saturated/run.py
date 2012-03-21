import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

from modules.direct_saturated.run import \
     (solve as direct_solve,
      extract_saturated_water_heights as extract_heights)

PARAMETERS = {
    'inverse-heights': { 'l1' : -1.0, 'wl1' : -1.0 }
    }

import modules.direct_saturated.run as direct

def adjust_cfg(flattened_cfg):
    direct.adjust_cfg(flattened_cfg)

mass_in_idx  = 0
mass_out_idx = 1

def solve(model):
    def ip_direct_saturated_heights(model, Ks):
        model.ks = Ks

        (_flag, t, z) = direct_solve(model)

        wl = extract_heights(z, model)

        return t, wl

    def lsq_ip_direct_saturated_heights(xdata, Ks):
        return ip_direct_saturated_heights(xdata, Ks)[1] # return only wl

    if model.exp_type == 'ish-f':
        # falling head test needs special adjustments
        model.include_acceleration = False
        # we ignore the r0 value just as omega in measured data (if present),
        # because it is of no validity for falling head test; instead we use
        # the (default) value in "r0-f" and compute the corresponding omega

        model.r0 = [model.r0_fall for wl0 in model.wl0]
        model.omega = (np.sqrt(model.g/model.r0)
                       * np.ones(np.shape(model.l0), float))

    # resolve the type of measured data
    if model.exp_type in ['ish', 'ish-sc', 'ish-f']:
        data_measured = model.wl1
        lsq_direct_fn = lsq_ip_direct_saturated_heights
        direct_fn     = ip_direct_saturated_heights
    elif model.exp_type in ['isc', 'isc-sc']:
        raise NotImplementedError('inverse::characteristics: load '
                                  ' characteristics data not implemented !')
        GC_measured   = data['GC']
        RM_measured   = data['RM']
        data_measured = np.concatenate((GC_measured, RM_measured))
        direct_fn    = ip_direct_saturated_characteristics

    else:
        raise ValueError('Unrecognized experiment type exp_type: %s'
                         % model.exp_type)
    xdata         = model

    model.ks = 1.1e-7
    Ks_init = model.ks # backup...
    # Solve inverse problem
    #    Ks_inv, cov_ks = curve_fit(lsq_direct_fn, xdata,
    #                           data_measured, p0 = Ks_init)
    Ks_inv, cov_ks = curve_fit(lsq_direct_fn, xdata,
                               data_measured, p0 = [Ks_init])

    t_inv, wl1_inv = direct_fn(xdata, Ks_inv)

    model.ks = Ks_init #...and restore

    # Print results
    for i in np.arange(len(data_measured)):
        print('Subexperiment ', i+1)
        print('    wl0:          % .6f' % model.wl0[i])
        print('    wl1_measured: % .6f    t_end_expected: % .2f' %
              (model.wl1[i],  model.duration[i]))
        print('    wl1_computed: % .6f    t_end_computed: % .2f' %
              (wl1_inv[i], t_inv[i]))
        print('    Error (%%):   % .2f                        % .2f' %
              ((wl1_inv[i] - model.wl1[i]) / model.wl1[i] * 100,
               (t_inv[i] - model.duration[i]) / model.duration[i] * 100))

    print('\nKs found: ', Ks_inv)

    return Ks_inv
