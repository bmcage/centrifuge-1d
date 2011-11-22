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
from direct_saturated import direct_saturated_problem, extract_saturated_characteristics

mass_in_idx  = 0
mass_out_idx = 1

[inifiles, outputdir, savecfgname] = \
    parse_centrifuge_input(sysargv)

def inverse_saturated(model, Ks):
    model.ks = Ks
    print(Ks)
    _flag, t, z = direct_saturated_problem(model)
    GC, RM = extract_saturated_characteristics(t, z, model)
    return np.concatenate((GC, RM))

def main():
    cfgmngr = ConfigManager.get_instance()
    model   = cfgmngr.read_configs(merge_output_p = True,
                                   preserve_sections_p = False,
                                   filenames = inifiles,
                                   defaults = [DEFAULT_PARAMETERS],
                                   saveconfigname = savecfgname)
    model.register_key('experiment', 'tspan', np.array([]))

    if not model.inverse_data_filename:
        raise ValueError('Data file for inverse problem not specified !')

    data = load_data(model.inverse_data_filename)
    t_measured  = data['t']
    GC_measured = data['GC']
    RM_measured = data['RM']
    model.tspan = t_measured
    data_measured = np.concatenate((GC_measured, RM_measured))
    data.close()

    Ks_init = model.ks
    measurements = np.concatenate((GC_measured, RM_measured))
    Ks_inv, cov_ks = curve_fit(inverse_saturated, model, data_measured, p0 = Ks_init)

    print('Ks initial: ', Ks_init)
    print('Ks found:   ', Ks_inv)
    print('Ks exact:   ', 2.4e-5)

if __name__ == "__main__":
    main()
