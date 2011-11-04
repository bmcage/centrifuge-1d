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
import matplotlib.pyplot as plt

from centrifugeparameters import CentrifugeParameters
from auxiliaryfunctions   import parse_centrifuge_input

syspath.append(join_path('.', 'odes', 'build', 'lib.linux-x86_64-3.2', 'scikits', 'odes', 'sundials'))
#/scikits/sundials_odes/build//
#print('syspath: ', syspath)
#import ida

from common_defs import ResFunction

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

def draw_graphs(fignum, t, wl_in, wl_out):
    plt.figure(fignum)
    plt.subplot(211)
    plt.plot(t, wl_in, 'b')
    plt.xlabel('Time')
    plt.ylabel('water in inflow chamber')

    plt.subplot(212)
    plt.plot(t, wl_out, 'k')
    plt.xlabel('Time')
    plt.ylabel('water in outflow chamber')
    plt.show()

class centrifuge_rhs(ResFunction):
    def evaluate(self, t, x, xdot, result, model):

        omega2g = ((model.omega_start  + (model.omega - model.omega_start)
                                  * (1 - np.exp(-model.omega_gamma*t)))**2 
                                  / model.g)
        L = model.l
        l = x[0]
        q_out = model.ks*omega2g * (l*(2*model.r0 - l)/L + 2*model.r0 + L)
        result[0] = xdot[0] + q_out
        result[1] = xdot[1] - q_out
        
        return 0
        
rhs = centrifuge_rhs()     

def main():
    cfgmngr = ConfigManager.get_instance()
    model   = cfgmngr.read_configs(merge_output_p = True,
                                   preserve_sections_p = False,
                                   filenames = inifiles,
                                   defaults = [DEFAULT_PARAMETERS],
                                   saveconfigname = savecfgname)
    model.register_key('experiment', 'tspan', \
                   arange(model.t_start, model.t_end, model.t_step))
    print(model.tspan)
    #try:
    import ida
    #except ImportError:
   #     print('ImportError: Could not load module ''ida''. Exiting...')
   #     return
    solver = ida.IDA()
    solver.set_options(resfn=rhs,
                       compute_initcond='yp0',
                       first_step=1e-18,
                       atol=1e-6,rtol=1e-6,
                       user_data=model)
    z0  = np.array([model.l0_in, 0])
    zp0 = np.zeros(z0.shape, float)
    z = solver.run_solver(model.tspan, z0, zp0)
    print(z)
    z = z[::10, :]
    tspan = tspan[::10, :]
    draw_graphs(1, tspan, z[:, 0], z[:, 1])

if __name__ == "__main__":
    main()
