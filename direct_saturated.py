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

from os.path import join as join_path
from sys import path as syspath
syspath.append(join_path('.', 'odes', 'build', 'lib.linux-x86_64-3.2'))

import scikits.odes.sundials.ida as ida
from scikits.odes.sundials.common_defs import ResFunction
from shared_functions import find_omega2g

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

class direct_saturated_rhs(ResFunction):
    def evaluate(self, t, x, xdot, result, model):

        omega2g = find_omega2g(t, model)
        L = model.l
        l = x[0]
        q_out = model.ks*omega2g * (l*(2*model.r0 - l)/L + 2*model.r0 + L)
        result[mass_in_idx] = xdot[0] + q_out
        result[mass_out_idx] = xdot[1] - q_out
        
        return 0

class direct_saturated_rhs_full(ResFunction):
    def evaluate(self, t, x, xdot, result, model):
        # x = [l_in, hL, l_out]
        omega2g = find_omega2g(t, model)
        L      = model.l
        l_in   = x[0]
        h0     = 1/2 * omega2g * l_in * (2*model.r0 - l_in)
        hL     = x[1]
        l_out  = x[2]
        qt     = model.ks * (omega2g/2. * (2*model.r0 + L) - (hL - h0) / L)

        #print('h0: ', h0, ' hL: ', hL, ', l_in: ', l_in, ', l_out: ', l_out, ', r0: ', model.r0)
        #print('x: ', x, 'xdot: ', xdot)
        print('x: ', x)
        #print('porosity*L: ', L*model.porosity, ', WM: ', model.water_volume,
        #      'l_in: ', l_in, ', l_out: ', l_out)
        result[0] = xdot[0] + qt
        result[1] = l_in + l_out + L*model.porosity - model.water_volume
        result[2] = xdot[2] - qt

        print('result: ', result)

        return 0

def total_water_volume(model):
    #TODO: we assume that we start with no outspelled water
    return (model.l0_in + model.l * model.porosity)

def extract_saturated_characteristics(t, z, model):
    GC, RM = characteristics(model.tspan, z[:, mass_in_idx],
                             z[:, mass_out_idx], model)[:2]

    return GC, RM

def extract_saturated_water_heights(z, model):
    return z[:, mass_in_idx]

rhs      = direct_saturated_rhs()
rhs_full = direct_saturated_rhs_full()

def solve_direct_saturated_problem(model):

    solver = ida.IDA(rhs,
                     compute_initcond='yp0',
                     first_step=1e-18,
                     atol=1e-6,rtol=1e-6,
                     user_data=model)
    z0  = np.array([model.l0_in, 0], float)
    zp0 = np.zeros(z0.shape, float)
    flag, t, z = solver.run_solver(model.tspan, z0, zp0)[:3]

    return flag, t, z

def solve_direct_saturated_problem_full(model):
    model.water_volume = total_water_volume(model)

    solver = ida.IDA(rhs_full,
                     compute_initcond='yp0',
                     first_step=1e-18,
                     atol=1e-6,rtol=1e-6,
                     algebraic_vars_idx=[1],
                     user_data=model)
    z0  = np.array([model.l0_in, 0, 0], float)
    zp0 = np.zeros(z0.shape, float)
    #solver.init_step(0, z0, zp0)
    flag, t, z = solver.run_solver(model.tspan, z0, zp0)[:3]

    return flag, t, z

def utilize_model(model):
    model.omega_start = model.omega_start / 60
    model.omega       = model.omega / 60

    model.register_key('experiment', 'tspan',
                       np.arange(model.t_start, model.t_end, model.t_step))
    model.register_key('experiment','water_volume', total_water_volume(model))

    return model

def run_direct(draw_graphs_p = False):
    from sys import argv as sysargv
    from auxiliaryfunctions import (parse_centrifuge_input,
                                    load_centrifuge_configs)
   
    [inifiles, outputdir, savecfgname] =  parse_centrifuge_input(sysargv)
    model = load_centrifuge_configs(inifiles, [utilize_model)
     
    _flag, t, z = solve_direct_saturated_problem(model)

    GC, RM = extract_saturated_characteristics(t, z, model)
    if draw_graphs_p:
        draw_graphs(1, t, z[:, 0], z[:, 1], GC, RM)
   
    return model, GC, RM

if __name__ == "__main__":
    run_direct(draw_graphs_p = True)
