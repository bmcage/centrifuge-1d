from __future__ import division

import numpy as np

from scikits.odes.sundials.common_defs import ResFunction
from modules.shared.functions import find_omega2g
from modules.shared.vangenuchten import h2Kh, u2Ku, dudh, dhdu, h2u, u2h
from modules.shared.solver import simulate_direct

class centrifuge_residual(ResFunction):

    def evaluate(self, t, w, wdot, result, model):
        # F(t,h,dh/dt) = porosity * du/dh * dh/dt
        #                + d[K(h)*(dh/dr - omega^2/g * r)]/dr

        last_idx      = model.last_idx
        saturated_idx = model.saturated_idx

        z    = w[:last_idx+1]

        # if positive pressure or negative saturation we want to refine the step
        if np.any(z[:saturated_idx] > 0.) or np.any(z[saturated_idx:] < 0.):
            return 1

        print('t ={:10.8f}'.format(t))

        inflow_p  = (wdot[model.mass_in_idx] > 0.)
        outflow_p = (wdot[model.s2_idx] == model.l0)

        s2 = w[model.s2_idx]
        s1 = w[model.s1_idx]
        ds = s2 - s1

        Ks = model.ks
        n  = model.n
        m  = model.m
        gamma = model.gamma
        rE  = model.re
        L   = model.l0
        fl2 = model.fl2
        r0  = rE - fl2 - L
        porosity = model.porosity


        h_saturated_idx = u2h(z[saturated_idx], n, m, gamma)

        z12  = (z[:-1] + z[1:]) / 2.
        z12[saturated_idx-1] = (z[saturated_idx-1] + h_saturated_idx) / 2.

        omega2g = find_omega2g(t, model)

        K12  = np.empty([model.inner_points + 1], dtype=float)
        h2Kh(z12[:saturated_idx], n, m, gamma, Ks, Kh=K12[:saturated_idx])
        u2Ku(z12[saturated_idx:],    m, gamma, Ks, Ku=K12[saturated_idx:])

        dy   = model.dy

        zdot =  wdot[:last_idx+1]

        dzdy12 = (z[1:] - z[:-1]) / model.dy
        dzdy12[saturated_idx -1] = ((h_saturated_idx - z[saturated_idx - 1])
                                    /model.dy[saturated_idx -1])

        dhdz12 = np.ones([model.inner_points + 1], dtype=float)
        dhdu(z12[saturated_idx:], n, m, gamma, dhdu=dhdz12[saturated_idx:])

        # TODO: repair q_first
        inflow_p = False
        if inflow_p:
            pass
        else:
            q_first = 0.

        q12 = -K12 * (dhdz12 * dzdy12 / ds - omega2g*(r0 + s1 + ds * model.y12))

        p = 3./2. + 1./m
        u_tmp = np.power(z[-3:], p)
        ds2dt  = -Ks*m*m/gamma/(n-1)/p/ds * (model.ldc1[-1] * u_tmp[0]
                                             - model.ldc2[-1] * u_tmp[1]
                                             + model.ldc3[-1] * u_tmp[2])


        dudz = np.ones([model.inner_points + 2], dtype=float)

        dudh(z12[:saturated_idx -1], n, m, gamma, dudh=dudz[:saturated_idx - 1])

        result[0] = (porosity * dudz[0] * zdot[0]
                             + 2 / dy[0] / ds * (q12[0] - q_first))

        result[1:last_idx] = (porosity * dudz[1:-1] * zdot[1:-1]
                                        + 2 / (dy[:-1] + dy[1:]) / ds
                                          * (q12[1:] - q12[:-1]))

        if outflow_p:
            pass
        else:
            result[last_idx]  = zdot[-1]
            result[model.mass_out_idx] = wdot[model.mass_out_idx]
            result[model.s2_idx]  = porosity*wdot[model.s2_idx] - ds2dt

        if inflow_p:
            pass
        else:
            result[model.mass_in_idx]  = wdot[model.mass_in_idx]
            result[model.s1_idx]  = wdot[model.s1_idx]


        result[model.pq_idx]  = wdot[model.pq_idx]

        return 0

    # TODO repair characteristics
def characteristics(t, u, mass_in, mass_out, s1, s2, model, chtype = 'all'):

    calc_gc = chtype in ['all', 'gc']
    calc_rm = chtype in ['all', 'rm']

    rE  = model.re
    L   = model.l0
    fl2 = model.fl2
    r0  = rE - fl2 - L

    porosity = model.porosity
    y  = model.y
    dy = model.dy
    l0_out = L + model.wt_out
    l_out  = L + model.wt_out - mass_out

    ds = s2 - s1

    if calc_rm:
        P = np.pi * model.d / 4
        omega2g = find_omega2g(t, model.omega, model)

    # Water mass
    wm_sat = ds/2  * (dy[0]* u[0] + dy[-1]*u[-1]
                      + np.sum((dy[:-1] + dy[1:])*u[1:-1]))

    WM    = model.density * (mass_in + porosity*wm_sat + mass_out)

    # GC is from the start of the sample (not from centr.axis)
    # sample = filter1 + soil + filter2
    r0_gc = model.fl1
    r0_rm = r0

     # Gravitational center
    r0 = r0_gc
    if calc_gc:
        gc_unsat = (porosity * 1/2 * model.density * ds
                    * ((r0 + s1)*dy[0]*u[0]
                       + (r0 + s2)*dy[-1]*u[-1]
                       + np.sum((dy[:-1]+dy[1:])
                                *(r0 + s1 + ds*y[1:-1])*u[1:-1])))
        gc_sat   = (1/2 * model.density
                    * (porosity * (np.power(r0 + s1, 2) - np.power(r0, 2))
                       + (np.power(r0, 2) - np.power(r0 - mass_in, 2))
                       + (np.power(r0 + l0_out, 2) - np.power(r0 + l_out, 2))))
        gc_left  = (gc_unsat + gc_sat) / WM
        gc_right = model.fl1 + model.l0 + model.fl2 - gc_left
        GC       = gc_right
    else:
        GC = None

    # Rotational momentum
    if calc_rm:
        r0 = r0_rm
        rm_unsat = (porosity * 1/4 * model.density * ds
                    * (np.power(r0 + s1, 2)*dy[0]*u[0]
                       + np.power(r0 + s2, 2)*dy[-1]*u[-1]
                       + np.sum((dy[:-1]+dy[1:]) * u[1:-1]
                                * np.power(r0 + s1 + ds*y[1:-1], 2))))
        rm_sat   = (1/6 * model.density
                    * (porosity * (np.power(r0 + s1, 3) - np.power(r0, 3))
                       + (np.power(r0, 3) - np.power(r0 - mass_in, 3))
                       + (np.power(r0 + l0_out, 3) - np.power(r0 + l_out, 3))))

        RM = omega2g * P * (rm_unsat + rm_sat)
    else:
        RM = None

    return GC, RM, WM

residual_fn = centrifuge_residual()

def solve(model):

    t   = np.empty([model.iterations+1, ], dtype=float)
    GC  = np.empty([model.iterations+1, ], dtype=float)
    WM  = np.empty([model.iterations+1, ], dtype=float)
    z   = np.empty([model.iterations+1, model.z_size], dtype=float)
    u   = np.empty([model.iterations+1, model.inner_points+2], dtype=float)
    h   = np.empty([model.iterations+1, model.inner_points+2], dtype=float)
    z0  = np.empty([model.z_size, ], float)

    saturated_idx = model.saturated_idx

    while model.next_iteration():
        i = model.iteration

        if i == 1:
            # initialize z0
            u0_first = h2u(model.h_init, model.n, model.m, model.gamma)
            u0_last  = 0.
            u0 = u0_last + (u0_first - u0_last) * (1-np.power(model.y, 2))

            u2h(u0[:saturated_idx], model.n, model.m, model.gamma,
                h=z0[:saturated_idx])
            z0[saturated_idx: model.last_idx+1] = \
              u0[saturated_idx: model.last_idx+1]

            s1 = 0.0
            s2 = model.l0 / 5.
            mass_in = mass_out = 0.0

            z0[model.s1_idx] = s1
            z0[model.s2_idx] = s2
            z0[model.mass_in_idx]  = mass_in
            z0[model.mass_out_idx] = mass_out

            # set values for t[0], z[0], u[0], GC[0]
            t[0] = 0.0
            z[0, :] = z0
            u[0, :] = u0

            GC[0] = characteristics(t[0], u0, mass_in, mass_out,
                                    s1, s2, model, chtype='gc')[0]
        else:
            z0 = z[i-1, :]

        (flag, t, z, i) = simulate_direct(model, residual_fn, z0)

        t[i] = t[i-1] + model.duration

        h2u(z[i, :saturated_idx], model.n, model.m, model.gamma,
            u=u[i, :saturated_idx])
        u[i, saturated_idx:] = z[i, saturated_idx: model.last_idx+1]

        u_min_limit = (-1e4 * model.gamma) ** (1-model.n)
        h[i, :] = z[i, :model.last_idx+1]
        h_view = h[i, saturated_idx:]

        h_view[h_view < u_min_limit] = u_min_limit
        u2h(h_view, model.n, model.m, model.gamma, h=h_view)

        s1 = z[i, model.s1_idx]
        s2 = z[i, model.s2_idx]
        #exit(0)
        mass_in  = z[i, model.mass_in_idx]
        mass_out = 0.0 # for GC and RM no expelled water is taken into account
        GC[i], _RM, WM[i] = characteristics(t_out, u[i, :], mass_in, mass_out,
                                            s1, s2, model, chtype='gc')

    if model.draw_graphs:
        from modules.shared.show import draw_graphs

        draw_graphs(t, y=model.y, s1=z[:, model.s1_idx], s2=z[:, model.s2_idx],
                    h=h, u=u, mass_out = z[:, model.mass_out_idx], GC=GC,
                    WM=WM)

    return (flag, t, z, GC)

def run(model):
    return solve(model)
