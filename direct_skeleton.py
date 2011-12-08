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

from sys import argv as sysargv, path as syspath
from os.path import join as join_path
from config import DEFAULT_PARAMETERS, ConfigManager
from numpy import arange
import numpy as np

from centrifugeparameters import CentrifugeParameters
from auxiliaryfunctions   import parse_centrifuge_input

syspath.append(join_path('.', 'odes', 'scikits', 'sundials_odes', 'build', 'lib.linux-x86_64-3.2'))
#/scikits/sundials_odes/build//
#print('syspath: ', syspath)
#import ida

from common_defs import ResFunction
from shared_functions import find_omega2g

[inifiles, outputdir, savecfgname] = \
    parse_centrifuge_input(sysargv)
    
print('savecfgname: ', savecfgname)

first_idx    = 0
last_idx     = -1
mass_in_idx  = -1
s1_idx       = -1
s2_idx       = -1
mass_out_idx = -1
pq_idx       = -1



class centrifuge_rhs(ResFunction):
    def evaluate(self, t, x, xdot, result, model):
        # porosity * du/dt = -Ks*(dh/dr - omega^2/g*r)
        # porosity * du/dh * dh/dt = -Ks*(dh/dr - omega^2/g*r)
        # F(t,h,dh/dt) = porosity * du/dh * dh/dt + Ks*(dh/dr - omega^2/g * r)
        
        omega2g = find_omega2g(t, model)

        result[0] = omega2g * x[0] * (2*model.r0 - result[0])
        return 0
        
rhs = centrifuge_rhs()     

def main():
    cfgmngr = ConfigManager.get_instance()
    model   = cfgmngr.read_configs(merge_output_p = True,
                                   preserve_sections_p = False,
                                   filenames = inifiles,
                                   defaults = [DEFAULT_PARAMETERS],
                                   saveconfigname = savecfgname)
    model.omega_start = model.omega_start / 60
    model.omega       = model.omega / 60
    
    #first_idx    = 0
    last_idx     = model.inner_points+1
    mass_in_idx  = model.inner_points+2
    s1_idx       = model.inner_points+3
    s2_idx       = model.inner_points+4
    mass_out_idx = model.inner_points+5
    pq_idx       = model.inner_points+6
    
    model.register_key('experiment', 'tspan', \
                   arange(model.t_start, model.t_end, model.t_step))
    if model.dtype == 1:
        x = np.linspace(0, 1, model.inner_points)
    else:
        raise NotImplemented('Currently only linear discretization is implemented')
    model.register_key('additional', 'x', x)
    model.register_key('additional', 'x12', (x[1:]+x[0:last_idx-1])/2.)
    
    
    try:
        import ida
    except ImportError:
        print('ImportError: Could not load module ''ida''. Exiting...')
        return
    solver = ida.IDA()
    solver.set_options(resfn=rhs,
                       compute_initcond='yp0',
                       first_step=1e-18,
                       atol=1e-6,rtol=1e-6,
                       algebraic_vars_idx=[4]
                       user_data=model)
    zp0 = np.zeros(model.z0.shape, float)
    solver.run_solver(model.tspan, model.z0, 

if __name__ == "__main__":
    main()
