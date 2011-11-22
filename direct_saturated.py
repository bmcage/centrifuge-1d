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

import numpy as np

from centrifugeparameters import CentrifugeParameters

from os.path import join as join_path
from sys import argv as sysargv, path as syspath
syspath.append(join_path('.', 'odes', 'build', 'lib.linux-x86_64-3.2'))

import scikits.odes.sundials.ida as ida
from scikits.odes.sundials.common_defs import ResFunction

mass_in_idx  = 0
mass_out_idx = 1

    
def draw_graphs(fignum, t, wl_in, wl_out, GC, RM):

    import matplotlib.pyplot as plt

    plt.figure(fignum, figsize=(12, 8))

    plt.subplot(221)
    plt.plot(t, wl_in, 'b')
    plt.xlabel('Time')
    plt.ylabel('Water in inflow chamber [cm]')

    plt.subplot(222)
    plt.plot(t, wl_out, 'k')
    plt.xlabel('Time')
    plt.ylabel('Water in outflow chamber [cm]')

    plt.subplot(223)
    plt.plot(t, wl_in, 'b')
    plt.xlabel('Time')
    plt.ylabel('Gravitational center [$cm$]')

    plt.subplot(224)
    plt.plot(t, wl_out, 'k')
    plt.xlabel('Time')
    plt.ylabel('Rotational momentum [$kg.m.s^{-2}]')

    plt.show()

def characteristics(t, mass_in, mass_out, model):
    porosity = model.porosity
    r0 = model.r0
    L  = model.l
    l0_out = L + model.l0_out
    l_out  = L + model.l0_out - mass_out

    omega2g = (np.power(model.omega_start  + (model.omega - model.omega_start)
                        * (1 - np.exp(-model.omega_gamma*t)), 2)
                / model.g)

    GC = np.empty(t.shape, float)
    RM = np.empty(t.shape, float)

    P = np.pi * model.d / 4

    WM = P * (mass_in + porosity * L + mass_out)

    for i in range(len(t)):
        # Gravitational center
        gc_sat   = (1/2 * model.density
                    * (porosity * (np.power(r0 + L, 2) - np.power(r0, 2))
                       + (np.power(r0, 2) - np.power(r0 - mass_in[i], 2))
                       + (np.power(r0 + l0_out, 2) - np.power(r0 + l_out[i], 2))))
        GC[i] =  P * gc_sat / WM[i]

        # Rotational momentum
        rm_sat   = (1/6 * model.density
                    * (porosity * (np.power(r0 + L, 3) - np.power(r0, 3))
                        + (np.power(r0, 3) - np.power(r0 - mass_in[i], 3))
                        + (np.power(r0 + l0_out, 3) - np.power(r0 + l_out[i], 3))))

        RM[i] = omega2g[i] * P * rm_sat

    return GC, RM, WM


class centrifuge_rhs(ResFunction):
    def evaluate(self, t, x, xdot, result, model):

        omega2g = ((model.omega_start  + (model.omega - model.omega_start)
                                  * (1 - np.exp(-model.omega_gamma*t)))**2 
                                  / model.g)
        L = model.l
        l = x[0]
        q_out = model.ks*omega2g * (l*(2*model.r0 - l)/L + 2*model.r0 + L)
        result[mass_in_idx] = xdot[0] + q_out
        result[mass_out_idx] = xdot[1] - q_out
        
        return 0


def extract_saturated_characteristics(t, z, model):
    GC, RM = characteristics(model.tspan, z[:, mass_in_idx],
                             z[:, mass_out_idx], model)[:2]

    return GC, RM

rhs = centrifuge_rhs()

def direct_saturated_problem(model):

    solver = ida.IDA(rhs,
                     compute_initcond='yp0',
                     first_step=1e-18,
                     atol=1e-6,rtol=1e-6,
                     user_data=model)
    z0  = np.array([model.l0_in, 0])
    zp0 = np.zeros(z0.shape, float)
    flag, t, z = solver.run_solver(model.tspan, z0, zp0)[:3]

    return flag, t, z

def run_direct(draw_graphs_p = False):
    from auxiliaryfunctions import parse_centrifuge_input
    from config import DEFAULT_PARAMETERS, ConfigManager
   
    [inifiles, outputdir, savecfgname] =  parse_centrifuge_input(sysargv)
    
    cfgmngr = ConfigManager.get_instance()
    model   = cfgmngr.read_configs(merge_output_p = True,
                                   preserve_sections_p = False,
                                   filenames = inifiles,
                                   defaults = [DEFAULT_PARAMETERS],
                                   saveconfigname = savecfgname)

     model.register_key('experiment', 'tspan',
                   np.arange(model.t_start, model.t_end, model.t_step))
     
     _flag, t, z = direct_saturated_problem(model.ks, model)

    GC, RM = extract_saturated_characteristics(t, z, model)
    if draw_graphs_p:
        draw_graphs(1, t, z[:, 0], z[:, 1], GC, RM)
   
    return GC, RM

if __name__ == "__main__":
    solve_direct(draw_graphs_p = True)
