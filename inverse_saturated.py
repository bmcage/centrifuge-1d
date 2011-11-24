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

from sys import argv as sysargv, path as syspath
from os.path import join as join_path
from config import DEFAULT_PARAMETERS, ConfigManager
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

from centrifugeparameters import CentrifugeParameters
from auxiliaryfunctions   import parse_centrifuge_input

syspath.append(join_path('.', 'odes', 'build', 'lib.linux-x86_64-3.2'))

from auxiliaryfunctions import load_data
from direct_saturated import direct_saturated_problem, extract_saturated_characteristics, extract_saturated_water_heights

mass_in_idx  = 0
mass_out_idx = 1

[inifiles, outputdir, savecfgname] = \
    parse_centrifuge_input(sysargv)

def inverse_saturated_characteristics(model, Ks):
    model.ks = Ks
    print(Ks)
    _flag, t, z = direct_saturated_problem(model)
    GC, RM = extract_saturated_characteristics(t, z, model)
    return np.concatenate((GC, RM))

def inverse_saturated_heights(xdata, Ks):
    model    = xdata[0]
    heights0 = xdata[1]
    model.ks = Ks
    print(Ks)
    #print(xdata)
    heights = np.empty(heights0.shape, float)
    for i in np.arange(len(heights0)):
        model.l0_in = heights0[i]
        _flag, t, z = direct_saturated_problem(model)
        h1 = extract_saturated_water_heights(z, model)
        print('%f -> %f' % (h1[0], h1[1]))
        heights[i] = h1[1]

    return heights

def main():
    cfgmngr = ConfigManager.get_instance()
    model   = cfgmngr.read_configs(merge_output_p = True,
                                   preserve_sections_p = False,
                                   filenames = inifiles,
                                   defaults = [DEFAULT_PARAMETERS],
                                   saveconfigname = savecfgname)
    model.register_key('experiment', 'tspan', np.array([]))
    model.omega_start = model.omega_start / 60
    model.omega       = model.omega / 60

    if not model.inverse_data_filename:
        raise ValueError('Data file for inverse problem not specified !')

    data = load_data(model.inverse_data_filename)

    t_measured    = data['t']
    model.tspan   = t_measured

    if model.inverse_data_type == 0:
        GC_measured   = data['GC']
        RM_measured   = data['RM']
        data_measured = np.concatenate((GC_measured, RM_measured))
        xdata         = model
        inverse_fn    = inverse_saturated_characteristics
    elif model.inverse_data_type == 1:
        heights_0     = data['heights0']
        heights_1     = data['heights1']
        data_measured = heights_1
        xdata         = (model, heights_0)
        inverse_fn    = inverse_saturated_heights
    else:
        raise ValueError('Unrecognized inverse_data_type: %d'
                         % model.inverse_data_type)
        
    data.close()

    #print(heights_0.dtype, heights_1.dtype, t_measured.dtype)
    Ks_init = model.ks
    Ks_inv, cov_ks = curve_fit(inverse_fn, xdata,
                               data_measured, p0 = Ks_init)

    print('Ks initial: ', Ks_init)
    print('Ks found:   ', Ks_inv)
    if model.inverse_data_type == 0:
        print('Ks exact:   ', 2.4e-5)
    elif model.inverse_data_type == 1:
        h1_inv = inverse_fn(xdata, Ks_inv)
        print('measured value - computed value')
        for i in np.arange(len(heights_0)):
            print('% f - % f' % (heights_1[i], h1_inv[i]))

if __name__ == "__main__":
    main()
