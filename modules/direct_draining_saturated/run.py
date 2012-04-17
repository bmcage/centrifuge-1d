import numpy as np

from scikits.odes.sundials.common_defs import ResFunction
from modules.shared.shared_functions import find_omega2g, h2Kh, dudh, h2u
from modules.shared.solver import simulate_direct

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

        omega2g = find_omega2g(t, model)

        s2 = model.l0
        s1 = 0
        ds = s2 - s1

        Ks = model.ks
        n  = model.n
        m  = model.m
        gamma = model.gamma
        r0 = model.r0
        L  = model.l0
        porosity = model.porosity

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
        q_last  = np.maximum(1e-12, -Kh_last * (dhdr_last/ds - omega2g*(r0 + L)))
        #print('ql:', q_last)
        #print('q_out', q_last)

        du_dh = dudh(h, n, m, gamma)
        result[first_idx] = (porosity * du_dh[0] * hdot[0]
                             + 2 / dy[0] / ds * (q12[0] - q_first))

        result[first_idx+1:last_idx] = (porosity * du_dh[1:-1] * hdot[1:-1]
                                        + 2 / (dy[:-1] + dy[1:]) / ds
                                          * (q12[1:] - q12[:-1]))
        result[last_idx]  = (porosity * du_dh[-1] * hdot[-1]
                             + 2 / dy[-1] / ds * (q_last - q12[-1]))

        result[model.mass_in_idx]  = zdot[model.mass_in_idx]
        result[model.mass_out_idx] = zdot[model.mass_out_idx]  - q_last
        result[model.s1_idx]  = zdot[model.s1_idx]
        result[model.s2_idx]  = zdot[model.s2_idx]
        result[model.pq_idx]  = zdot[model.pq_idx]

        return 0

def characteristics(t, u, mass_in, mass_out, s1, s2, model, chtype = 'all'):

    calc_gc = chtype in ['all', 'gc']
    calc_rm = chtype in ['all', 'rm']

    porosity = model.porosity
    y  = model.y
    dy = model.dy
    L  = model.l0
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
    r0_rm = model.r0

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
    z   = np.empty([model.iterations+1, model.z_size], dtype=float)
    u   = np.empty([model.iterations+1, model.inner_points+2], dtype=float)
    z0  = np.empty([model.z_size, ], float)

    model.init_iteration()

    while model.next_iteration():
        i = model.iteration

        if i == 1:
             # initialize: z0
             # set values for t[0], z[0], u[0], GC[0]
             z0[model.first_idx:model.last_idx+1] = model.h_init
             s1 = 0.0
             s2 = model.l0
             mass_in = mass_out = 0.0

             z0[model.s1_idx] = s1
             z0[model.s2_idx] = s2
             z0[model.mass_in_idx]  = mass_in
             z0[model.mass_out_idx] = mass_out

             t[0] = 0.0
             z[0, :] = z0
             u[0, :] = h2u(z0[model.first_idx: model.last_idx+1],
                           model.n, model.m, model.gamma)

             GC[0] = characteristics(t[0], u[0, :], mass_in, mass_out,
                                     s1, s2, model, chtype='gc')[0]
        else:
            z0 = z[i-1, :]

        (flag, t_out, z[i, :]) = simulate_direct(model, residual_fn, z0)

        t[i] = t[i-1] + model.duration

        u[i, :] = h2u(z[i, model.first_idx: model.last_idx+1],
                      model.n, model.m, model.gamma)
        s1 = z[i, model.s1_idx]
        s2 = z[i, model.s2_idx]

        mass_in  = z[i, model.mass_in_idx]
        mass_out = 0.0 # for GC and RM no expelled water is taken into account
        GC[i] = characteristics(t_out, u[i, :], mass_in, mass_out,
                                          s1, s2, model, chtype='gc')[0]

    if model.draw_graphs:
        from modules.shared.show import draw_graphs

        h = z[:, model.first_idx:model.last_idx+1]
        GC[:] = 0.5

        draw_graphs(t, y=model.y, h=h, u=u, mass_out=z[:, model.mass_out_idx],
                    GC=GC, s1=z[:, model.s1_idx], s2=z[:, model.s2_idx])

    return (flag, t, z, GC)
