from sys import path as syspath
from os.path import join as join_path
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

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

def ip_direct_saturated_heights(xdata, Ks):
    model       = xdata[0]
    heights0    = xdata[1]
    return_time = xdata[2]
    model.ks = Ks

    heights = np.empty(heights0.shape, float)
    t_end   = np.empty(heights0.shape, float)
    for i in np.arange(len(heights0)):
        model.l0_in = heights0[i]
        _flag, t, z = solve_direct_saturated_problem(model)

        h1 = extract_saturated_water_heights(z, model)
        heights[i] = h1[1]
        t_end[i]   = t[1] # t = [start_time, end_time]

        if model.debugging:
            print('t = ', t, ', Ks =  ', Ks, ', h0 = ', model.l0_in,
                  ', h_found = ', heights[i], ', r0 = ', model.r0,
                  'omega = ', model.omega, ' flag = ', _flag)
    if return_time:
        return heights, t_end
    else:
        return heights

def solve_inverse_saturated(model, measured_data_filename):
    data = load_measured_data(measured_data_filename)

    t_duration     = np.asarray(data.duration, float)
    t_deceleration = model.deceleration_duration
    model.t_start  = 0.
    model.t_end    = t_duration
    if model.include_acceleration:
        model.tspan    = np.array([0, t_duration + model.deceleration_duration], float)
    else:
         model.tspan    = np.array([0, t_duration], float)
    model.omega    = data.omega
    model.l        = data.length
    model.r0       = data.r0

    if model.data_type == 0:
        raise NotImplemented('inverse::characteristics: load characteristics data not implemented !')
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
        xdata         = (model, heights_0, return_time)

        direct_fn = ip_direct_saturated_heights

        if model.data_type == 2:# falling head test
            model.include_acceleration = False
            model.omega = np.sqrt(model.g / model.r0)
    else:
        raise ValueError('Unrecognized data_type: %d'
                         % model.data_type)

    Ks_init = model.ks
     
    Ks_inv, cov_ks = curve_fit(direct_fn, xdata,
                               data_measured, p0 = Ks_init)

    return_time   = True # return also times for estimating error
    xdata         = (model, heights_0, return_time)
    h1_inv, t_inv = direct_fn(xdata, Ks_inv)

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

def run_inverse_saturated():
    from sys import argv as sysargv
    from auxiliaryfunctions import (parse_centrifuge_input,
                                    load_centrifuge_configs)

    inifiles = parse_centrifuge_input(sysargv)[0]
    model = load_centrifuge_configs(inifiles,
                                    [verify_inverse_data, utilize_direct_model])

    data_filenames = model.inverse_data_filename
    Ks_init = model.ks
    Ks_inv  = np.empty([len(data_filenames), ],  float)
    
    for i in range(len(data_filenames)):
        model.ks  = Ks_init
        Ks_inv[i] = solve_inverse_saturated(model, data_filenames[i])

    Ks_inv_disp = Ks_inv / 100
    print('\n\nKs_inv computed: ')
    for Ks_print in Ks_inv_disp:
        print(Ks_print)

    return Ks_inv


if __name__ == "__main__":
    run_inverse_saturated()
