import numpy as np

from os.path import join as join_path
from sys import path as syspath
syspath.append(join_path('.', 'odes', 'build', 'lib.linux-x86_64-3.2'))

import scikits.odes.sundials.ida as ida
from scikits.odes.sundials.common_defs import ResFunction
from shared_functions import find_omega2g

mass_in_idx  = 0
mass_out_idx = 1

    
def draw_graphs(fignum, t, wl_in, wl_out, GC = None, RM = None):

    import matplotlib.pyplot as plt

    if (not GC is None) or (not RM is None):
        rows = 2
        figheight = 8 
    else:
        rows = 1
        figheight = 4

    plt.figure(fignum, figsize=(12, figheight))

    plt.subplot(rows,2,1)
    plt.plot(t, wl_in, 'b')
    plt.xlabel('Time')
    plt.ylabel('Water in inflow chamber [cm]')

    plt.subplot(rows,2,2)
    plt.plot(t, wl_out, 'k')
    plt.xlabel('Time')
    plt.ylabel('Water in outflow chamber [cm]')

    if not GC is None:
        plt.subplot(rows,2,3)
        plt.plot(t, wl_in, 'b')
        plt.xlabel('Time')
        plt.ylabel('Gravitational center [$cm$]')

    if not RM is None:
        plt.subplot(rows,2,4)
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

    GC = np.empty(t.shape, float)
    RM = np.empty(t.shape, float)

    P = np.pi * model.d / 4

    WM = P * (mass_in + porosity * L + mass_out)

    for i in range(len(t)):
        omega2g = find_omega2g(t, model)
        
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

        RM[i] = omega2g * P * rm_sat

    return GC, RM, WM

class direct_saturated_rhs(ResFunction):
    def evaluate(self, t, x, xdot, result, model):

        omega2g = find_omega2g(t, model)
        L = model.l
        l = x[0]
        q_out = model.ks*omega2g/2. * (l*(2*model.r0 - l)/L + 2*model.r0 + L)
        result[mass_in_idx] = xdot[0] + q_out
        result[mass_out_idx] = xdot[1] - q_out
        
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

direct_saturated_rhs_fn = direct_saturated_rhs()

def solve_direct_saturated_problem(model):

    solver = ida.IDA(direct_saturated_rhs_fn,
                     compute_initcond='yp0',
                     first_step=1e-18,
                     atol=1e-6,rtol=1e-6,
                     user_data=model)
    z0  = np.array([model.l0_in, 0], float)
    zp0 = np.zeros(z0.shape, float)
    flag, t, z = solver.run_solver(model.tspan, z0, zp0)[:3]

    return flag, t, z

def utilize_model(model):
    model.register_key('experiment','water_volume', total_water_volume(model))

def run_direct(draw_graphs_p = False):
    from sys import argv as sysargv
    from auxiliaryfunctions import (parse_centrifuge_input,
                                    load_centrifuge_configs)
   
    [inifiles, outputdir, savecfgname] =  parse_centrifuge_input(sysargv)
    model = load_centrifuge_configs(inifiles, [utilize_model])
     
    _flag, t, z = solve_direct_saturated_problem(model)

    if model.data_type == 0:
        GC, RM = extract_saturated_characteristics(t, z, model)
        if draw_graphs_p:
            draw_graphs(1, t, z[:, 0], z[:, 1], GC, RM)
        return model, GC, RM
    elif model.data_type == 1:
        h1 = extract_saturated_water_heights(z, model)
        print('t:  ', t)
        print('h1: ', h1)
        if draw_graphs_p:
            draw_graphs(1, t, z[:, 0], z[:, 1])
        
        return h1
    else:
        raise ValueError('direct_saturated::run_direct Unknown data type: ', data_type)

if __name__ == "__main__":
    run_direct(draw_graphs_p = True)
