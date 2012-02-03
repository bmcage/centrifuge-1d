from sys import path as syspath
from os.path import join as join_path
import numpy as np
import matplotlib.pyplot as plt

syspath.append(join_path('.', 'odes', 'build', 'lib.linux-x86_64-3.2'))

import scikits.odes.sundials.ida as ida
from scikits.odes.sundials.common_defs import ResFunction
from shared_functions import find_omega2g, h2Kh, dudh, h2u, characteristics

first_idx    = 0
last_idx     = -1
mass_in_idx  = -1
s1_idx       = -1
s2_idx       = -1
mass_out_idx = -1
pq_idx       = -1
z_size       = -1

def draw_graphs(fignum, t, y, z, model):
    h = z[:, first_idx:last_idx+1]
    u = h2u(h, model.n, model.m, model.gamma)
    s1 = z[:, s1_idx]
    s2 = z[:, s2_idx]
    ds = s2 - s1
    mass_in  = z[:, mass_in_idx]
    mass_out = z[:, mass_out_idx]
    x = np.empty([len(t), len(y)], float)

    GC, RM, MW = characteristics(t, u, mass_in, mass_out, s1, s2, model)

    legend_data = []
    for i in range(len(t)):
        x[i, :] = z[i, s1_idx] + y * ds[i]
        legend_data.append(str('t = % d' % t[i]))

    plt.figure(fignum, figsize=(16, 8.5))

    plt.subplot(321)
    plt.plot(x.transpose(), h.transpose(), '.')
    plt.xlabel('Rotational axis distance ''r'' [$cm$]')
    plt.ylabel('Piezometric head ''h'' [$cm$]')

    plt.subplot(322)
    plt.plot(x.transpose(), u.transpose(), '.')
    plt.xlabel('Rotational axis distance ''r'' [$cm$]')
    plt.ylabel('Relative saturation ''u''')
    plt.legend(legend_data, bbox_to_anchor=(1., 1.), loc=2, borderaxespad=0.1)

    plt.subplot(323)
    plt.plot(t, z[:, mass_out_idx], '.')
    plt.xlabel('Time [$s$]')
    plt.ylabel('Outspelled water [$cm^3$]')

    plt.subplot(325)
    plt.plot(t, GC.transpose(), '.')
    plt.xlabel('Time [$s$]')
    plt.ylabel('Gravitational center [$cm$]')

    plt.subplot(326)
    plt.plot(t, RM.transpose(), '.')
    plt.xlabel('Time [$s$]')
    plt.ylabel('Rotational momentum [$kg.m.s^{-1}$]')
    
    plt.show()

class centrifuge_residual(ResFunction):
         
    def evaluate(self, t, z, zdot, result, model):
        # F(t,h,dh/dt) = porosity * du/dh * dh/dt
        #                + d[K(h)*(dh/dr - omega^2/g * r)]/dr

        h    =  z[first_idx:last_idx+1]

        if np.any(h > 0): # positive pressure - we want to refine the step
            return 1
        
        hdot =  zdot[first_idx:last_idx+1]
        h12  = (h[1:] + h[:-1]) / 2

        omega2g = find_omega2g(t, model)
        
        s2 = model.l
        s1 = 0
        ds = s2 - s1

        Ks = model.ks
        n  = model.n
        m  = model.m
        gamma = model.gamma
        r0 = model.r0
        L  = model.l

        Kh12 = h2Kh(h12, n, m, gamma, Ks)
        Kh_last =  h2Kh(h[-1], n, m, gamma, Ks)

        dy   = model.dy
        
        dhdr12 = (h[1:] - h[:-1]) / model.dy
        dhdr_last = (model.ldc1[-1] * h[-3]
                     - model.ldc2[-1] * h[-2]
                     + model.ldc3[-1] * h[-1])

        q_first = 0.
        q12 = -Kh12 * (dhdr12 / ds - omega2g*(r0 + s1 + ds * model.y12))
        q_last  = -Kh_last * np.minimum(0., dhdr_last/ds - omega2g*(r0 + L))

        du_dh = dudh(h, n, m, gamma, Ks)
        result[first_idx] = (model.porosity * du_dh[0] * hdot[0]
                             + 2 / dy[0] / ds * (q12[0] - q_first))
                             
        result[first_idx+1:last_idx] = (model.porosity * du_dh[1:-1] * hdot[1:-1]
                             + 2 / (dy[:-1] + dy[1:]) / ds * (q12[1:] - q12[:-1]))
        result[last_idx]  = (model.porosity * du_dh[-1] * hdot[-1]
                             + 2 / dy[-1] / ds * (q_last - q12[-1]))

        result[mass_in_idx]  = zdot[mass_in_idx]
        result[mass_out_idx] = zdot[mass_out_idx]  - q_last
        result[s1_idx]  = zdot[s1_idx]
        result[s2_idx]  = zdot[s2_idx]
        result[pq_idx]  = zdot[pq_idx]
        
        return 0

        
residual_fn = centrifuge_residual()

def utilize_model(model):
    from shared_functions import (lagrangean_derivative_coefs,  dudh, h2Kh,
                                  h2u, characteristics, find_omega2g)
    
    model.register_key('additional', 'm', 1-1/model.n)

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

def solve_direct_draining_saturated_problem(model):
    global z_size
    z0  = np.zeros([z_size, ], float)
    z0[first_idx:last_idx+1] = -0.15
    z0[s2_idx] = model.l
    zp0 = np.zeros([z_size, ], float)
    
    solver = ida.IDA(residual_fn,
                     #compute_initcond='yp0',
                     first_step=1e-20,
                     atol=1e-1,rtol=1e-2,
                     #algebraic_vars_idx=[4],
                     user_data=model)

    print('tsp',model.tspan)
    flag, t, z = solver.solve(model.tspan, z0, zp0)[:3]

    return flag, t, z

def run_direct_draining_saturated():
    from sys import argv as sysargv
    from auxiliaryfunctions import (parse_centrifuge_input,
                                    load_centrifuge_configs)

    inifiles = parse_centrifuge_input(sysargv)[0]
    model    = load_centrifuge_configs(inifiles,
                                       [utilize_model])

    global first_idx, last_idx, mass_in_idx, s1_idx, s2_idx
    global mass_out_idx, pq_idx, z_size
    #first_idx    = 0
    last_idx     = model.inner_points+1
    mass_in_idx  = model.inner_points+2
    s1_idx       = model.inner_points+3
    s2_idx       = model.inner_points+4
    mass_out_idx = model.inner_points+5
    pq_idx       = model.inner_points+6

    z_size       = model.inner_points+7 # total length of 'z' array

    
    flag, t, z = solve_direct_draining_saturated_problem(model)
    draw_graphs(1, t, model.y, z, model)

if __name__ == "__main__":
    run_direct_draining_saturated()
