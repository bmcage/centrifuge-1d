# -*- coding: utf-8 -*-
#
# Copyright (C) 2011-12  Pavol Kišon
# Copyright (C) 2009-12  Benny Malengier
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful, 
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
#

from sys import path as syspath
from os.path import join as join_path
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

syspath.append(join_path('.', 'odes', 'build', 'lib.linux-x86_64-3.2'))

from auxiliaryfunctions import load_measured_data
from direct_saturated import (solve_direct_saturated_problem,
                              solve_direct_saturated_problem_full,
                              extract_saturated_characteristics,
                              extract_saturated_water_heights,
                              utilize_model as utilize_direct_model)

mass_in_idx  = 0
mass_out_idx = 1



def inverse_saturated_characteristics(model, Ks):
    model.ks = Ks
    if model.debugging:
        print(Ks)
    _flag, t, z = solve_direct_saturated_problem(model)
    GC, RM = extract_saturated_characteristics(t, z, model)
    return np.concatenate((GC, RM))

def inverse_saturated_heights(xdata, Ks):
    model    = xdata[0]
    heights0 = xdata[1]
    model.ks = Ks

    heights = np.empty(heights0.shape, float)
    for i in np.arange(len(heights0)):
        model.l0_in = heights0[i]
        _flag, t, z = solve_direct_saturated_problem(model)

        h1 = extract_saturated_water_heights(z, model)
        heights[i] = h1[1]

        if model.debugging:
            print('t = ', t, ', Ks =  ', Ks, ', h0 = ', model.l0_in,
                  ', h_found = ', heights[i], ', r0 = ', model.r0,
                  'omega = ', model.omega, ' flag = ', _flag)

    return heights

def inverse_saturated_heights_full(xdata, Ks):
    model    = xdata[0]
    heights0 = xdata[1]
    model.ks = Ks

    heights = np.empty(heights0.shape, float)

    for i in np.arange(len(heights0)):
        model.l0_in = heights0[i]
        _flag, t, z = solve_direct_saturated_problem_full(model)
        h1 = extract_saturated_water_heights(z, model)
        print('%f -> %f' % (h1[0], h1[1]))
        heights[i] = h1[1]

    return heights

def solve_inverse_saturated(model, measured_data_filename):
    data = load_measured_data(measured_data_filename)

    t_duration    = np.asarray(data.duration, float)
    model.tspan   = np.array([0, t_duration], float)
    model.omega   = data.omega
    model.l       = data.length

    if model.inverse_data_type == 0:
        raise NotImplemented('inverse::characteristics: load characteristics data not implemented !')
        GC_measured   = data['GC']
        RM_measured   = data['RM']
        data_measured = np.concatenate((GC_measured, RM_measured))
        xdata         = model
        inverse_fn    = inverse_saturated_characteristics
    elif model.inverse_data_type == 1 or model.inverse_data_type == 2:
        heights_0     = np.asarray(data.h0, float)
        heights_1     = np.asarray(data.h1, float)
        data_measured = heights_1
        xdata         = (model, heights_0)
        if model.inverse_data_type == 1:
            inverse_fn = inverse_saturated_heights
        else:
            raise NotImplemented('inverse: inverse_data_type =2 (full problem) not implemented')
            inverse_fn = inverse_saturated_heights_full
    else:
        raise ValueError('Unrecognized inverse_data_type: %d'
                         % model.inverse_data_type)

    Ks_init = model.ks
     
    Ks_inv, cov_ks = curve_fit(inverse_fn, xdata,
                               data_measured, p0 = Ks_init)

    h1_inv = inverse_fn(xdata, Ks_inv)
    if model.debugging:
        print('measured value - computed value')
    for i in np.arange(len(heights_0)):
        print('% f - % f' % (heights_1[i], h1_inv[i]))

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
                         % model.inverse_data_type)

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
    print('Ks_inv: ')
    for Ks_print in Ks_inv_disp:
        print(Ks_print)

    return Ks_inv


if __name__ == "__main__":
    run_inverse_saturated()
