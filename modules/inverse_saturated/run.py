from sys import path as syspath
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

from modules.direct_saturated.run import \
     (solve as direct_solve,
      extract_saturated_water_heights as extract_heights)
from config import merge_cfgs

syspath.append('/'.join(['.', 'odes', 'build', 'lib.linux-x86_64-3.2']))

INVERSE_SATURATED_ADDITIONAL_PARAMETERS = {
    'inverse-heights': { 'l1' : -1.0, 'wl1' : -1.0 }
    }

import modules.direct_saturated.run as direct

def base_cfg():
    return merge_cfgs(direct.base_cfg(), INVERSE_SATURATED_ADDITIONAL_PARAMETERS)

def adjust_cfg(flattened_cfg):
    direct.adjust_cfg(flattened_cfg)

mass_in_idx  = 0
mass_out_idx = 1

def check_cfg(flattened_cfg):
    #TODO
    return True

# def ip_direct_saturated_characteristics(model, Ks):
#     model.ks = Ks
#     if model.debugging:
#         print(Ks)
#     _flag, t, z = solve_direct_saturated_problem(model)
#     GC, RM = extract_saturated_characteristics(t, z, model)
#     return np.concatenate((GC, RM))




# add to check_fn:
#if model.l0 < 0:
# if np.isscalar(data.omega):
#            if data.omega < 0:
 # if np.isscalar(data.r0):
 #            if data.r0 < 0:
 #                xdata_r0 = np.asarray([model.r0], float)
 #            else:
 #                xdata_r0 = np.asarray([data.r0], float)
 #        else:
 #            xdata_r0 = np.asarray(data.r0, float)
 # if np.isscalar(data.duration):
 #        if data.t_duration < 0:
 
    # else:
    #     xdata_t_end = np.asarray(data.duration, float)

    # xdata_t_start  = np.zeros(np.shape(xdata_t_end), float)

    # xdata_t_span = np.zeros([np.alen(xdata_t_end), 2], float)

    # if model.include_acceleration:
    #     xdata_t_span[:, 1] = xdata_t_end + model.deceleration_duration
    # else:
    #     xdata_t_span[:, 1] = xdata_t_end

def solve(model):
    def ip_direct_saturated_heights(model, Ks):
        #print('omg:',model.omega)
        #result_heights = np.empty(model.wl0.shape, float)
        #result_t_end   = np.empty(result_heights, float)
        model.ks = Ks
        #print(Ks)
        _flag, t, z = direct_solve(model)
        #print(z)
        wl = extract_heights(z, model)
        #result_heights[i] = h1[1]
        #result_t_end[i]   = t[1] # t = [start_time, end_time]

        #return  result_t_end, result_heights
        return t, wl

    def lsq_ip_direct_saturated_heights(xdata, Ks):
        return ip_direct_saturated_heights(xdata, Ks)[1] # return only wl

    if model.exp_type == 'ish-f':
        # falling head test needs special adjustments
        model.include_acceleration = False
        # we ignore the r0 value just as omega in measured data (if present),
        # because it is of no validity for falling head test; instead we use 
        # the (default) value in "r0-f" and compute the corresponding omega

        #print('1. data.omega: ', data.omega, 'data.r0: ', data.r0,
        #      'model.r0', model.r0)
        model.r0 = model.r0_f
        model.omega = (np.sqrt(model.g/model.r0)
                       * np.ones(np.shape(model.l0), float))
        #print('2. data.omega: ', data.omega, 'data.r0: ', data.r0,
        #      'model.r0', model.r0, 'xdataa_omega: ', xdata_omega, 
        #      'xdata_r0: ', xdata_r0)


        #print(model.exp_type, model.wl1, model.ks)
    # resolve the type of measured data
    if model.exp_type in ['ish', 'ish-sc']:
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
        raise ValueError('Unrecognized data_type: %d'
                         % model.data_type)
    xdata         = model

    model.ks = 1.1e-8
    Ks_init = model.ks # backup...
    # Solve inverse problem
    #    Ks_inv, cov_ks = curve_fit(lsq_direct_fn, xdata,
    #                           data_measured, p0 = Ks_init)
    Ks_inv, cov_ks = curve_fit(lsq_direct_fn, xdata,
                               data_measured, p0 = [Ks_init])

    t_inv, wl1_inv = direct_fn(xdata, Ks_inv)
    #print('tspan: ',model.tspan)
    model.ks = Ks_init #...and restore

    # Print results
    for i in np.arange(len(data_measured)):
        print('Subexperiment ', i+1)
        print('    wl0:          % .6f' % model.wl0[i])
        print('    wl1_measured: % .6f    t_end_expected: % .2f' %
              (model.wl1[i],  model.tspan[1]))
        print('    wl1_computed: % .6f    t_end_computed: % .2f' %
              (wl1_inv[1], t_inv[1]))
        print('    Error (%%):   % .2f                        % .2f' %
              ((wl1_inv[1] - model.wl1[i]) / model.wl1[i] * 100,
               (t_inv[1] - model.tspan[1]) / model.tspan[1] * 100))

        print('\nKs found: ', Ks_inv)

    return Ks_inv

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

def unify_data(model):
    pass

# def run_inverse_saturated():
#     from sys import argv as sysargv
#     from auxiliaryfunctions import (parse_centrifuge_input,
#                                     load_centrifuge_configs)

#     inifiles = parse_centrifuge_input(sysargv)[0]
#     imodel = load_centrifuge_configs(inifiles,
#                                      [verify_inverse_data, utilize_direct_model])

#     data = [load_measured_data(filename) for filename in model.inverse_data_filename]

#     # save original values of model and replace them with read data
#     backup_include_acceleration = imodel.include_acceleration
#     backup_t_start  = imodel.t_start
#     backup_t_end    = imodel.t_end
#     backup_t_span   = imodel.tspan
#     backup_omega    = imodel.omega
#     backup_l        = imodel.l
#     backup_r0       = imodel.r0
#     Ks_init         = imodel.ks

#     Ks_inv  = np.empty([len(data_filenames), ],  float)

#     for i in range(len(data)):
#         Ks_inv[i] = solve_inverse_saturated(imodel, data_filenames[i])

#         # Restore original values
#         model.include_acceleration = backup_include_acceleration
#         model.t_start = backup_t_start
#         model.t_end   = backup_t_end
#         model.t_span  = backup_t_span
#         model.omega   = backup_omega
#         model.l       = backup_l
#         model.r0      = backup_r0
#         model.ks      = Ks_init

#     Ks_inv_disp = Ks_inv / 100
#     print('\n\nKs_inv computed [m/s]: ')
#     for Ks_print in Ks_inv_disp:
#         print(Ks_print)

#     return Ks_inv


if __name__ == "__main__":
    #run_inverse_saturated()
    pass
