import numpy as np

from sys import path as syspath
syspath.append('/'.join(['.', 'odes', 'build', 'lib.linux-x86_64-3.2']))

import scikits.odes.sundials.ida as ida
from scikits.odes.sundials.common_defs import ResFunction
#import cProfile

from shared_functions import find_omega2g, h2Kh, dudh, h2u
from common import set_model_attributes

PARAMETERS = {'fluid': ['density'],
     'discretization': ['inner_points', 'dtype'],
         'centrifuge': ['wt_out']}

# wt_out - distance of the output chamber (that catches outspelled water)
#          from the end of the sample (tube)

CFG_ADDITIONAL_PARAMETERS = {}

def base_cfg():
    from modules.base.run import base_cfg as raw_cfg
    from config import merge_cfgs

    return merge_cfgs(raw_cfg(), CFG_ADDITIONAL_PARAMETERS)

def adjust_cfg(flattened_cfg):
    from  modules.base.run import adjust_cfg as base_adjust_cfg
    base_adjust_cfg(flattened_cfg)

    for param in ['r0_fall', 'wt_out', 'inner_points', 'dtype', 'density',
                  'h_init', 'n', 'gamma', 'draw_graphs']:
        if not param in flattened_cfg:
            raise ValueError('CFG:check: Value not initialized: %s' % param)

    from shared_functions import (lagrangean_derivative_coefs,  dudh, h2Kh,
                                  h2u, characteristics, find_omega2g)

    gamma = flattened_cfg['gamma']
    if len(gamma) > 1:
        print('Parameter gamma is a list of length > 1. Only constant or '
              'list of length 1 is supported. Exiting...')
        exit(1)
    else:
        flattened_cfg['gamma'] = gamma[0]

    n = flattened_cfg['n']
    if type(n) == list:
        if len(n) > 1:
            print('Parameter n is a list of length > 1. Only constant or '
                  'list of length 1 is supported. Exiting...')
            exit(1)
        else:
            flattened_cfg['n'] = n[0]

    flattened_cfg['m'] = 1-1/flattened_cfg['n']

    inner_points = flattened_cfg['inner_points']

    flattened_cfg['first_idx']    = 0
    flattened_cfg['last_idx']     = inner_points+1
    flattened_cfg['mass_in_idx']  = inner_points+2
    flattened_cfg['s1_idx']       = inner_points+3
    flattened_cfg['s2_idx']       = inner_points+4
    flattened_cfg['mass_out_idx'] = inner_points+5
    flattened_cfg['pq_idx']       = inner_points+6

    # total length of 'z' array (discretization points + s1,s2,mass_in,...)
    flattened_cfg['z_size']       = inner_points+7 

    if flattened_cfg['dtype'] == 1:
        y = np.linspace(0, 1, inner_points + 2)
    else:
        raise NotImplementedError('For now only linear discretization is'
                                  ' implemented')

    flattened_cfg['y']   = y
    flattened_cfg['y12'] = (y[1:]+y[:-1])/2.
    
    dy = y[1:]-y[:-1]
    flattened_cfg['dy']  = dy
    
    ldc1, ldc2, ldc3 = lagrangean_derivative_coefs(dy)
    flattened_cfg['ldc1']  = ldc1
    flattened_cfg['ldc2']  = ldc2
    flattened_cfg['ldc3']  = ldc3

    if ('wl0' in flattened_cfg and
        ((type(flattened_cfg['wl0']) == list and any(flattened_cfg['wl0']))
         or flattened_cfg['wl0'] != 0.)):
        raise ValueError("InputArgError: 'wl0' is specified and is not zero.\n"
                         "For drainage saturated non-zero 'wl0' is not allowed.")

    flattened_cfg['_l0'] = -1.0

def draw_graphs(fignum, t, y, z, model):
    import matplotlib.pyplot as plt
    
    h = z[:, model.first_idx:model.last_idx+1]
    u = h2u(h, model.n, model.m, model.gamma)
    s1 = z[:, model.s1_idx]
    s2 = z[:, model.s2_idx]
    ds = s2 - s1
    mass_in  = z[:, model.mass_in_idx]
    mass_out = 0.0 # for GC no expelled water is taken into account
    x = np.empty([len(t), len(y)], float)

    GC, _RM, _MW = characteristics(t, u, mass_in, mass_out, s1, s2, model,
                                   chtype='gc')

    legend_data = []
    for i in range(len(t)):
        x[i, :] = z[i, model.s1_idx] + y * ds[i]
        legend_data.append('t =%7d' % t[i])

    plt.figure(fignum, figsize=(16, 8.5))

    plt.subplots_adjust(wspace=0.15, left=0.06, right=0.85)

    plt.subplot(321)
    plt.plot(x.transpose(), h.transpose(), '.')
    plt.xlabel('Rotational axis distance ''r'' [$cm$]')
    plt.ylabel('Piezometric head ''h'' [$cm$]')

    plt.subplot(322)
    plt.plot(x.transpose(), u.transpose(), '.')
    plt.xlabel('Rotational axis distance ''r'' [$cm$]')
    plt.ylabel('Relative saturation ''u''')
    plt.legend(legend_data, bbox_to_anchor=(1.02, 1.), loc=2, borderaxespad=0.0)

    plt.subplot(323)
    plt.plot(t, z[:, model.mass_out_idx], '.')
    plt.xlabel('Time [$s$]')
    plt.ylabel('Outspelled water [$cm^3$]')

    plt.subplot(325)
    plt.plot(t, GC.transpose(), '.')
    plt.xlabel('Time [$s$]')
    plt.ylabel('Gravitational center [$cm$]')

    # plt.subplot(326)
    # plt.plot(t, RM.transpose(), '.')
    # plt.xlabel('Time [$s$]')
    # plt.ylabel('Rotational momentum [$kg.m.s^{-1}$]')
    
    plt.show()

class centrifuge_residual(ResFunction):
         
    def evaluate(self, t, z, zdot, result, model):
        # F(t,h,dh/dt) = porosity * du/dh * dh/dt
        #                + d[K(h)*(dh/dr - omega^2/g * r)]/dr

        first_idx = model.first_idx
        last_idx  = model.last_idx

        h    =  z[first_idx:last_idx+1]

        if np.any(h > 0): # positive pressure - we want to refine the step
            return 1
        #print('h: ', h)
        
        hdot =  zdot[first_idx:last_idx+1]
        h12  = (h[1:] + h[:-1]) / 2

        omega2g = find_omega2g(t, model._omega, model)

        s2 = model._l0
        s1 = 0
        ds = s2 - s1

        Ks = model.ks
        n  = model.n
        m  = model.m
        gamma = model.gamma
        r0 = model._r0
        L  = model._l0
        porosity = model._porosity

        Kh12 = h2Kh(h12, n, m, gamma, Ks)
        Kh_last =  h2Kh(h[-1], n, m, gamma, Ks)

        #Kh12 = cProfile.run('h2Kh(h12, n, m, gamma, Ks)', 'h2kh1')
        #Kh_last =  cProfile.run('h2Kh(h[-1], n, m, gamma, Ks)', 'h2kh2')

        dy   = model.dy
        
        dhdr12 = (h[1:] - h[:-1]) / model.dy
        dhdr_last = (model.ldc1[-1] * h[-3]
                     - model.ldc2[-1] * h[-2]
                     + model.ldc3[-1] * h[-1])

        q_first = 0.
        q12 = -Kh12 * (dhdr12 / ds - omega2g*(r0 + s1 + ds * model.y12))
        q_last  = -Kh_last * np.minimum(0., dhdr_last/ds - omega2g*(r0 + L))
        #print('ql:', q_last)
        #print('q_out', q_last)

        du_dh = dudh(h, n, m, gamma, Ks)
        result[first_idx] = (porosity * du_dh[0] * hdot[0]
                             + 2 / dy[0] / ds * (q12[0] - q_first))

        result[first_idx+1:last_idx] = (porosity * du_dh[1:-1] * hdot[1:-1]
                             + 2 / (dy[:-1] + dy[1:]) / ds * (q12[1:] - q12[:-1]))
        result[last_idx]  = (porosity * du_dh[-1] * hdot[-1]
                             + 2 / dy[-1] / ds * (q_last - q12[-1]))

        result[model.mass_in_idx]  = zdot[model.mass_in_idx]
        result[model.mass_out_idx] = zdot[model.mass_out_idx]  - q_last
        result[model.s1_idx]  = zdot[model.s1_idx]
        result[model.s2_idx]  = zdot[model.s2_idx]
        result[model.pq_idx]  = zdot[model.pq_idx]
        
        return 0

def characteristics(t, u, mass_in, mass_out, s1, s2, model, chtype='all'):
    calc_gc = (chtype == 'all' or chtype == 'gc')
    calc_rm = (chtype == 'all' or chtype == 'rm')

    porosity = model._porosity
    y  = model.y
    dy = model.dy
    L  = model._l0
    l0_out = L + model.wt_out
    l_out  = L + model.wt_out - mass_out

    ds = s2 - s1

    if calc_gc: 
        GC = np.empty(t.shape, float)
    else:
        GC = np.asarray([], dtype=float)

    if calc_rm: 
        RM = np.empty(t.shape, float)
        P = np.pi * model.d / 4
        
        omega2g = (np.power(model.omega_start
                            + (model.omega - model.omega_start)
                            * (1 - np.exp(-model.omega_gamma*t)), 2)
                / model.g)

    else:
        RM = np.asarray([], dtype=float)


    # Water mass
    wm_sat = ds/2  * (dy[0]* u[:, 0] + dy[-1]*u[:, -1]
                      + np.sum((dy[:-1] + dy[1:])*u[:, 1:-1], 1))

    WM = model.density * (mass_in + porosity*wm_sat + mass_out)

    # GC is from the start of the sample (not from centr.axis)
    # sample = filter1 + soil + filter2
    r0_gc = model.fl1
    r0_rm = model._r0
    for i in range(len(t)):
        # Gravitational center
        r0 = r0_gc
        if calc_gc:
            gc_unsat = (porosity * 1/2 * model.density * ds[i]
                        * ((r0 + s1[i])*dy[0]*u[i, 0]
                            + (r0 + s2[i])*dy[-1]*u[i, -1]
                            + np.sum((dy[:-1]+dy[1:])
                                     *(r0 + s1[i] + ds[i]*y[1:-1])*u[i, 1:-1])))
            gc_sat   = (1/2 * model.density
                        * (porosity * (np.power(r0 + s1[i], 2) - np.power(r0, 2))
                           + (np.power(r0, 2) - np.power(r0 - mass_in[i], 2))
                           + (np.power(r0 + l0_out, 2) - np.power(r0 + l_out, 2))))
            gc_left  = (gc_unsat + gc_sat) / WM[i]
            print(model.fl1, model._l0, model.fl2, gc_left)
            gc_right = model.fl1 + model._l0 + model.fl2 - gc_left
            GC[i]    = gc_right
            #print('GC: unsat, sat, gc/wm', gc_unsat, gc_sat, GC[i])

        # Rotational momentum
        if calc_rm:
            r0 = r0_rm
            rm_unsat = (porosity * 1/4 * model.density * ds[i]
            * (np.power(r0 + s1[i], 2)*dy[0]*u[i, 0]
               + np.power(r0 + s2[i], 2)*dy[-1]*u[i, -1]
               + np.sum((dy[:-1]+dy[1:]) * u[i, 1:-1]
                        * np.power(r0 + s1[i] + ds[i]*y[1:-1], 2))))
            rm_sat   = (1/6 * model.density
                        * (porosity * (np.power(r0 + s1[i], 3) - np.power(r0, 3))
                           + (np.power(r0, 3) - np.power(r0 - mass_in[i], 3))
                           + (np.power(r0 + l0_out, 3) - np.power(r0 + l_out[i], 3))))

            RM[i] = omega2g[i] * P * (rm_unsat + rm_sat)

    return GC, RM, WM
        
residual_fn = centrifuge_residual()

def solve(model):

    def run_solve(model):
        z0  = np.zeros([model.z_size, ], float)
        z0[model.first_idx:model.last_idx+1] = model.h_init # initial pressure h
        z0[model.s2_idx] = model._l0
        zp0 = np.zeros([model.z_size, ], float)

        #print('z0:', z0)
    
        solver = ida.IDA(residual_fn,
                     #compute_initcond='yp0',
                     first_step_size=1e-20,
                     atol=1e-1,rtol=1e-2,
                     #algebraic_vars_idx=[4],
                     #linsolver='band', uband=1, lband=1,
                     user_data=model)

        flag, t, z = solver.solve(model.tspan, z0, zp0)[:3]

        return flag, t, z

    t = np.empty([len(model.duration)+1, ], dtype=float)
    z = np.empty([len(t), model.z_size], dtype=float) # 2 columns: wl_in, wl_out

    if not (model.ks1 or model.fl1):
        model.ks1 = -1.0
        model.fl1 =  0.0
    if not (model.ks2 or model.fl2):
        model.ks2 = -1.0
        model.fl2 =  0.0

   

    solver = ida.IDA(residual_fn,
                     #compute_initcond='yp0',
                     first_step_size=1e-20,
                     atol=1e-1,rtol=1e-2,
                     max_step_size=840.,
                     max_steps=8000,
                     #algebraic_vars_idx=[4],
                     linsolver='band', uband=1, lband=2,
                     user_data=model)

    tspans = np.cumsum(model.duration)

    attributes_list =  ['l0', 'r0', 'omega', 'porosity']

    for i in range(len(model.duration)):

        # if model.ks1 and model.fl1:
        #     model._ks1 = model.ks1[i]
        #     model._fl1 = model.fl1[i]
        # if model.ks2 and model.fl2:
        #     model._ks2 = model.ks2[i]
        #     model._fl2 = model.fl2[i]

        #print(model.l0, model.r0, model.wl0)
        #print('l0: ', model.l0, 'r0: ', model.r0)        
        #print('ks2: ', model.ks2)
        #print('rf:', model.r0_fall)

        set_model_attributes(model, i, attributes_list)

        #flag, tacc, zacc  = run_solve(model)
        if i == 0:
             z0  = np.zeros([model.z_size, ], float)
             z0[model.first_idx:model.last_idx+1] = model.h_init # initial pressure h
             z0[model.s2_idx] = model._l0
             zp0 = np.zeros([model.z_size, ], float)
             t[0] = 0.0
             z[0, :] = z0

             solver._init_step(0.0, z0, zp0)
             
        flag, t_out = solver.step(tspans[i], z[i+1, :])
        t[i+1] = t_out

        #t[i]    = tacc[1]
        #z[i, :] = zacc[1, :]

        # falling head
        if hasattr(model, 'fh_duration') and model.fh_duration[i] > 0.:
            acc = model.include_acceleration
            model.tspan  = np.asarray([0.0, model.fh_duration[i]])
            model._r0    = model.r0_fall
            model._omega = model.omega_fall[i]

            flag, tfh, zfh = run_solve(model)
            model.include_acceleration = acc

    #print('t,z: ', t, z)

    if model.draw_graphs:
        draw_graphs(1, t, model.y, z, model)

    return (flag, t, z)

def run_direct_draining_saturated():
    from sys import argv as sysargv
    from auxiliaryfunctions import (parse_centrifuge_input,
                                    load_centrifuge_configs)

    inifiles = parse_centrifuge_input(sysargv)[0]
    model    = load_centrifuge_configs(inifiles,
                                       [utilize_model])

    
    
    flag, t, z = solve_direct_draining_saturated_problem(model)
    draw_graphs(1, t, model.y, z, model)

if __name__ == "__main__":
    run_direct_draining_saturated()
