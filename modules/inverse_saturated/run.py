from sys import path as syspath
from os.path import join as join_path
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

import modules.direct_saturated.direct_saturated

syspath.append(join_path('.', 'odes', 'build', 'lib.linux-x86_64-3.2'))

from auxiliaryfunctions import load_measured_data
from direct_saturated import (solve_direct_saturated_problem,
                              extract_saturated_characteristics,
                              extract_saturated_water_heights,
                              utilize_model as utilize_direct_model)

mass_in_idx  = 0
mass_out_idx = 1



def ip_direct_saturated_characteristics(model, Ks):
    model.ks = Ks
    if model.debugging:
        print(Ks)
    _flag, t, z = solve_direct_saturated_problem(model)
    GC, RM = extract_saturated_characteristics(t, z, model)
    return np.concatenate((GC, RM))

def lsq_ip_direct_saturated_heights(xdata, Ks):
    return ip_direct_saturated_heights(xdata, Ks)[1] # return only heights

def ip_direct_saturated_heights(xdata, Ks):
    model       = xdata[0]
    heights0    = xdata[1]
    t_start     = xdata[2]
    t_end       = xdata[3]
    t_span      = xdata[4]
    omega       = xdata[5]
    L           = xdata[6]
    r0          = xdata[7]
    model.ks = Ks

    result_heights = np.empty(heights0.shape, float)
    result_t_end   = np.empty(heights0.shape, float)

    for i in np.arange(len(heights0)):
        model.l0_in = heights0[i]
        model.t_start = t_start[i]
        model.t_end   = t_end[i]
        model.tspan   = t_span[i, :]
        model.omega   = omega[i]
        model.l       = L[i]
        model.r0      = r0[i]
        
        _flag, t, z = solve_direct_saturated_problem(model)

        h1 = extract_saturated_water_heights(z, model)
        result_heights[i] = h1[1]
        result_t_end[i]   = t[1] # t = [start_time, end_time]

        if model.debugging:
            print('t = ', t, ', Ks =  ', Ks, ', h0 = ', model.l0_in,
                  ', h_found = ', heights[i], ', r0 = ', model.r0,
                  'omega = ', model.omega, ' flag = ', _flag)

    return  result_t_end, result_heights


def solve_inverse_saturated(model, measured_data_filename):
    



    if np.isscalar(data.length):
        if data.length < 0:
            xdata_length = np.asarray([model.length], float)
        else:
            xdata_length = np.asarray([data.length], float)
    else:
        xdata_length = np.asarray(data.length, float)

    if model.data_type == 2:# falling head test needs special adjustments
        model.include_acceleration = False
        # we ignore the r0 value just as omega in measurements data (if present),
        # because it is of no validity (for falling head test); instead we use 
        # the values supplied in the original config file (which is in the
        # 'model' supplied

        #print('1. data.omega: ', data.omega, 'data.r0: ', data.r0,
        #      'model.r0', model.r0)
    
        xdata_omega = (np.sqrt(model.g/model.r0)
                       * np.ones(np.shape(xdata_length), float))
        xdata_r0    = model.r0 * np.ones(np.shape(xdata_length), float)
        #print('2. data.omega: ', data.omega, 'data.r0: ', data.r0,
        #      'model.r0', model.r0, 'xdataa_omega: ', xdata_omega, 
        #      'xdata_r0: ', xdata_r0)
    
    else:
        if np.isscalar(data.omega):
            if data.omega < 0:
                xdata_omega = np.asarray([model.omega], float)
            else:
                xdata_omega = np.asarray([data.omega], float)
        else:
            xdata_omega = np.asarray(data.omega, float)

        if np.isscalar(data.r0):
            if data.r0 < 0:
                xdata_r0 = np.asarray([model.r0], float)
            else:
                xdata_r0 = np.asarray([data.r0], float)
        else:
            xdata_r0 = np.asarray(data.r0, float)

    if np.isscalar(data.duration):
        if data.t_duration < 0:
            xdata_t_end = np.asarray([model.t_end], float)
        else:
            xdata_t_end = np.asarray([data.duration], float)
    else:
        xdata_t_end = np.asarray(data.duration, float)

    xdata_t_start  = np.zeros(np.shape(xdata_t_end), float)

    xdata_t_span = np.zeros([np.alen(xdata_t_end), 2], float)
    if model.include_acceleration:
        xdata_t_span[:, 1] = xdata_t_end + model.deceleration_duration
    else:
        xdata_t_span[:, 1] = xdata_t_end

    # resolve the type of measured data
    if model.data_type == 0:
        raise NotImplementedError('inverse::characteristics: load '
                                  ' characteristics data not implemented !')
        GC_measured   = data['GC']
        RM_measured   = data['RM']
        data_measured = np.concatenate((GC_measured, RM_measured))
        xdata         = model
        direct_fn    = ip_direct_saturated_characteristics
    elif model.data_type == 1 or model.data_type == 2:
        heights_0     = np.asarray(data.h0, float)
        heights_1     = np.asarray(data.h1, float)
        data_measured = heights_1
        return_time   = False
        xdata         = (model, heights_0, xdata_t_start, xdata_t_end,
                         xdata_t_span, xdata_omega, xdata_length, xdata_r0,
                         return_time)

        lsq_direct_fn = lsq_ip_direct_saturated_heights
        direct_fn     = ip_direct_saturated_heights
    else:
        raise ValueError('Unrecognized data_type: %d'
                         % model.data_type)

    # Solve inverse problem
    Ks_inv, cov_ks = curve_fit(lsq_direct_fn, xdata,
                               data_measured, p0 = Ks_init)

    t_inv, h1_inv = direct_fn(xdata, Ks_inv)

    # Print results
    for i in np.arange(len(heights_0)):
        print('Subexperiment ', i+1)
        print('    h0:          % .6f' % heights_0[i])
        print('    h1_measured: % .6f    t_end_expected: % .2f' %
              (heights_1[i],  model.tspan[1]))
        print('    h1_computed: % .6f    t_end_computed: % .2f' %
              (h1_inv[i], t_inv[i]))
        print('    Error (%%):   % .2f                        % .2f' %
              ((h1_inv[i] - heights_1[i]) / heights_1[i] * 100,
               (t_inv[i] - model.tspan[1]) / model.t_end * 100))

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

def run_inverse_saturated():
    from sys import argv as sysargv
    from auxiliaryfunctions import (parse_centrifuge_input,
                                    load_centrifuge_configs)

    inifiles = parse_centrifuge_input(sysargv)[0]
    imodel = load_centrifuge_configs(inifiles,
                                     [verify_inverse_data, utilize_direct_model])

    data = [load_measured_data(filename) for filename in model.inverse_data_filename]

    # save original values of model and replace them with read data
    backup_include_acceleration = imodel.include_acceleration
    backup_t_start  = imodel.t_start
    backup_t_end    = imodel.t_end
    backup_t_span   = imodel.tspan
    backup_omega    = imodel.omega
    backup_l        = imodel.l
    backup_r0       = imodel.r0
    Ks_init         = imodel.ks

    Ks_inv  = np.empty([len(data_filenames), ],  float)

    for i in range(len(data)):
        Ks_inv[i] = solve_inverse_saturated(imodel, data_filenames[i])

        # Restore original values
        model.include_acceleration = backup_include_acceleration
        model.t_start = backup_t_start
        model.t_end   = backup_t_end
        model.t_span  = backup_t_span
        model.omega   = backup_omega
        model.l       = backup_l
        model.r0      = backup_r0
        model.ks      = Ks_init

    Ks_inv_disp = Ks_inv / 100
    print('\n\nKs_inv computed [m/s]: ')
    for Ks_print in Ks_inv_disp:
        print(Ks_print)

    return Ks_inv


if __name__ == "__main__":
    run_inverse_saturated()
