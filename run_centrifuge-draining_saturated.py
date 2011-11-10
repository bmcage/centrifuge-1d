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

syspath.append(join_path('.', 'odes', 'build', 'lib.linux-x86_64-3.2'))

from scikits.odes.sundials.common_defs import ResFunction

[inifiles, outputdir, savecfgname] = parse_centrifuge_input(sysargv)
    
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

def lagrangean_derivative_coefs(dx):
    """
    Returns the coeficients for the Lagrangeand derivative of the differences
    array 'dx'. The first point has a right derivative, last point a left
    derivative and central difference is for the mid-points.
    """
    ldc1 = np.concatenate(([-(2*dx[0]+dx[1])/(dx[0]*(dx[0]+dx[1]))],
                          -dx[1:]/(dx[:-1]*(dx[:-1]+dx[1:])),
                          [dx[-1]/(dx[-2]*(dx[-2]+dx[-1]))]))
    ldc2 = np.concatenate(([dx[0]+dx[1]/(dx[1]*dx[0])],
                          (dx[1:] - dx[:-1])/dx[:-1]/dx[1:],
                          [(dx[-1]+dx[-2])/(dx[-2]*dx[-1])]))
    ldc3 = np.concatenate(([-dx[0]/(dx[1]*(dx[1]+dx[0]))],
                           dx[:-1]/(dx[1:]*(dx[:-1]+dx[1:])),
                           [(2*dx[-1]+dx[-2])/(dx[-1]*(dx[-2]+dx[-1]))]))

    return ldc1, ldc2, ldc3

def h2Kh(h, n, m, gamma, Ks):
    xkh  = np.power( gamma * h, n-1.)
    xkh2 = np.power(1 + gamma * h * xkh, m/2.)
    Kh   = Ks/xkh2 * np.power(1-xkh/np.power(xkh2, 2), 2)

    return Kh

def dudh(h, n, m, gamma, Ks):
    xkh  = np.power(gamma*h, n-1)
    xkh2 = np.power(1 + gamma*h * xkh, m+1)
    dudh = - gamma*(n-1) * xkh / xkh2

    return dudh
        
class centrifuge_rhs(ResFunction):
         
    def evaluate(self, t, z, zdot, result, model):
        # F(t,h,dh/dt) = porosity * du/dh * dh/dt
        #                + d[K(h)*(dh/dr - omega^2/g * r)]/dr
        
        omega2g = ((model.omega_start  + (model.omega - model.omega_start)
                                  * (1 - np.exp(-model.omega_gamma*t)))**2 
                                  / model.g)
        s2 = model.l
        s1 = 0
        ds = s2 - s1

        Ks = model.ks
        n  = model.n
        m  = model.m
        gamma = model.gamma
        r0 = model.r0
        L  = model.l
        
        h    =  z[first_idx:last_idx+1]
        hdot =  zdot[first_idx:last_idx+1]
        h12  = (h[1:] + h[:-1]) / 2

        Kh12 = h2Kh(h12, n, m, gamma, Ks)
        Kh_last =  h2Kh(h[-1], n, m, gamma, Ks)

        dy   = model.dy
        
        dhdr12 = (h[1:] - h[:-1]) / model.dy
        dhdr_last = (model.ldc1[-1] * h[-3]
                     - model.ldc2[-1] * h[-2]
                     + model.ldc3[-1] * h[-1])

        q_first = 0.
        q12 = Kh12 * (dhdr12 / ds - omega2g*(r0 + s1 + ds * model.y12))
        q_last  = Kh_last * np.minimum(0., dhdr_last/ds - omega2g*(r0 + L))

        du_dh = dudh(h, n, m, gamma, Ks)
        result[first_idx] = (model.porosity * du_dh[0] * hdot[0]
                             + 2 / dy[0] / ds * (q12[0] - q_first))
                             
        result[first_idx+1:last_idx] = (model.porosity * du_dh[1:-1] * hdot[1:-1]
                             + 2 / (dy[:-1] + dy[1:]) / ds * (q12[1:] - q12[:-1]))
        result[last_idx]  = (model.porosity * du_dh[-1] * hdot[-1]
                             + 2 / dy[-1] / ds * (q_last - q12[0]))

        result[mass_in_idx]  = 0.
        result[mass_out_idx] = -q_last
        result[s1_idx]  = 0.
        result[s2_idx]  = 0.
        result[pq_idx]  = 0.

        # print('dudh = ', du_dh)
        print('result[1:5] = ', result[:5])
        print('result[-10:] = ', result[-10:])
        
        return 0
        
rhs = centrifuge_rhs()     

def main():
    cfgmngr = ConfigManager.get_instance()
    model   = cfgmngr.read_configs(merge_output_p = True,
                                   preserve_sections_p = False,
                                   filenames = inifiles,
                                   defaults = [DEFAULT_PARAMETERS],
                                   saveconfigname = savecfgname)
    
    global first_idx, last_idx, mass_in_idx, s1_idx, s2_idx, mass_out_idx, pq_idx
    #first_idx    = 0
    last_idx     = model.inner_points+1
    mass_in_idx  = model.inner_points+2
    s1_idx       = model.inner_points+3
    s2_idx       = model.inner_points+4
    mass_out_idx = model.inner_points+5
    pq_idx       = model.inner_points+6

    z_size       = model.inner_points+7 # total length of 'z' array

    model.register_key('additional', 'm', 1-1/model.n)
    model.register_key('experiment', 'tspan',
                   arange(model.t_start, model.t_end, model.t_step))
    if model.dtype == 1:
        y = np.linspace(0, 1, model.inner_points+2)
    else:
        raise NotImplemented('Currently only linear discretization is implemented')
    model.register_key('additional', 'y', y)
    model.register_key('additional', 'y12', (y[1:]+y[:-1])/2.)
    dy = y[1:]-y[:-1]
    model.register_key('additional', 'dy', dy)
    ldc1, ldc2, ldc3 = lagrangean_derivative_coefs(dy)
    model.register_key('additional', 'ldc1', ldc1)
    model.register_key('additional', 'ldc2', ldc2)
    model.register_key('additional', 'ldc3', ldc3)

    try:
        from scikits.odes.sundials import ida
    except ImportError:
        print('ImportError: Could not load module ''ida''. Exiting...')
        return

    z0  = np.zeros([z_size, ], float)
    z0[first_idx:last_idx+1] = -100.
    zp0 = np.zeros([z_size, ], float)
    
    solver = ida.IDA(rhs,
                     #compute_initcond='yp0',
                     first_step=1e-18,
                     atol=1e-6,rtol=1e-2,
                     #algebraic_vars_idx=[4],
                     user_data=model)

    print(model.tspan)
    #    solver.run_solver(model.tspan, z0, zp0)
    z_ic0  = np.zeros([z_size, ], float)
    zp_ic0 = np.zeros([z_size, ], float)
    
    solver.init_step(model.tspan[0], z0, zp0, z_ic0, zp_ic0)
    solver.step(model.tspan[1],z_ic0, zp_ic0)

    print(z_ic0, zp_ic0)

if __name__ == "__main__":
    main()
