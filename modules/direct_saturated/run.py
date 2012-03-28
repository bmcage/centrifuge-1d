import numpy as np

import scikits.odes.sundials.ida as ida
from scikits.odes.sundials.common_defs import ResFunction
from shared_functions import find_omega2g

mass_in_idx  = 0
mass_out_idx = 1

def adjust_cfg(flattened_cfg):
    from  base import adjust_cfg as base_adjust_cfg
    base_adjust_cfg(flattened_cfg)

    flattened_cfg['_r0']    = -1.0
    flattened_cfg['_omega'] = -1.0
    flattened_cfg['_l0']    = -1.0
    flattened_cfg['_ks1']   = -1.0
    flattened_cfg['_ks2']   = -1.0
    flattened_cfg['_fl1']   =  0.0
    flattened_cfg['_fl2']   =  0.0

    flattened_cfg['fh_duration'] = np.asarray([], dtype=float)

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

class direct_saturated_rhs(ResFunction):
    def evaluate(self, t, x, xdot, result, model):

        omega2g = find_omega2g(t, model._omega, model)
        #print('o2g', omega2g)
        L = model._l0
        wl = x[0]
        rS = model._r0 - model._fl1 - wl
        rE = model._r0 + L + model._fl2

        qt = (omega2g/2. / (model._fl1/model._ks1 + L/model.ks
                            + model._fl2/model._ks2)
              * (rE*rE - rS*rS))

        #print('L', L, x, model._r0, model._fl1, model._ks1, model._fl2, model._ks2)
        #print('qt: ', qt, x, xdot, result)
        #print('t: ', t, 'x: ', x)
        result[mass_in_idx]  = xdot[0] + qt
        result[mass_out_idx] = xdot[1] - qt

        return 0

def extract_saturated_water_heights(z, model):
    return z[:, mass_in_idx]

direct_saturated_rhs_fn = direct_saturated_rhs()

def solve(model):

    def run_solve(model, tspan, z0):
        solver = ida.IDA(direct_saturated_rhs_fn,
                         compute_initcond='yp0',
                         first_step_size=1e-18,
                         atol=1e-6,rtol=1e-6,
                         user_data=model)
        zp0 = np.zeros(z0.shape, float)
        flag, t, z = solver.solve(tspan, z0, zp0)[:3]

        return flag, t, z

    t = np.empty([len(model.duration), ], dtype=float)
    z = np.empty([len(t), 2], dtype=float) # 2 columns: wl_in, wl_out

    if not (model.ks1 or model.fl1):
        model._ks1 = -1.0
        model._fl1 =  0.0
    if not (model.ks2 or model.fl2):
        model._ks2 = -1.0
        model._fl2 =  0.0

    tspan = np.asarray([0.0, 0.0], dtype=float)
    z0    = np.asarray([0.0, 0.0], dtype=float)

    for i in range(len(model.duration)):

        if model.ks1 and model.fl1:
            model._ks1 = model.ks1[i]
            model._fl1 = model.fl1[i]
        if model.ks2 and model.fl2:
            model._ks2 = model.ks2[i]
            model._fl2 = model.fl2[i]
        model._l0    = model.l0[i]
        model._r0    = model.r0[i]

        # acceleration
        tspan[1]  = model.duration[i]
        z0[0]     = model.wl0[i]
        model._omega = model.omega[i]

        flag, tacc, zacc  = run_solve(model, tspan, z0)

        t[i]    = tacc[1]
        z[i, :] = zacc[1, :]

        # falling head
        if model.fh_duration:
            tspan[1]     = model.fh_duration[i]
            z0[0]        = zacc[1, 0]
            model._r0    = model.r0_fall
            model._omega = model.omega_fall

            flag, tfh, zfh = run_solve(model, tspan, z0)

    return (flag, t, z)

    # TODO: finish with set options
