import numpy as np

from scikits.odes.sundials.common_defs import ResFunction
from modules.shared.shared_functions import find_omega2g
from modules.shared.vangenuchten import h2Kh, dudh, h2u
from modules.shared.characteristics import water_mass, calc_gc, calc_rm
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

        #Kh12 = cProfile.run('h2Kh(h12, n, m, gamma, Ks)', 'h2kh1')
        #Kh_last =  cProfile.run('h2Kh(h[-1], n, m, gamma, Ks)', 'h2kh2')

        dy   = model.dy

        dhdr12 = (h[1:] - h[:-1]) / model.dy

        q_first = 0.
        q12 = -Kh12 * (dhdr12 / ds - omega2g*(r0 + s1 + ds * model.y12))

        #print('ql:', q_last)
        #print('q_out', q_last)

        du_dh = dudh(h, n, m, gamma)
        result[first_idx] = (porosity * du_dh[0] * hdot[0]
                             + 2 / dy[0] / ds * (q12[0] - q_first))

        result[first_idx+1:last_idx] = (porosity * du_dh[1:-1] * hdot[1:-1]
                                        + 2 / (dy[:-1] + dy[1:]) / ds
                                          * (q12[1:] - q12[:-1]))

        if model.rb_type == 0:
            q_last = 0.
            result[last_idx]  = (porosity * du_dh[-1] * hdot[-1]
                                 + 2 / dy[-1] / ds * (q_last - q12[-1]))
        elif model.rb_type == 1:
            dhdr_last = (model.ldc1[-1] * h[-3]
                         - model.ldc2[-1] * h[-2]
                         + model.ldc3[-1] * h[-1])
            Kh_last =  h2Kh(h[-1], n, m, gamma, Ks)
            q_last  = np.maximum(1e-12,
                                 -Kh_last * (dhdr_last/ds - omega2g*(r0 + L)))
            result[last_idx]  = (porosity * du_dh[-1] * hdot[-1]
                                 + 2 / dy[-1] / ds * (q_last - q12[-1]))
        else:
            dhdr_last = (model.ldc1[-1] * h[-3]
                         - model.ldc2[-1] * h[-2]
                         + model.ldc3[-1] * h[-1])
            Kh_last =  h2Kh(h[-1], n, m, gamma, Ks)
            q_last  = np.maximum(1e-12,
                                 -Kh_last * (dhdr_last/ds - omega2g*(r0 + L)))
            result[last_idx]  = hdot[-1]

        result[model.mass_in_idx]  = zdot[model.mass_in_idx]
        result[model.mass_out_idx] = zdot[model.mass_out_idx]  - q_last
        result[model.s1_idx]  = zdot[model.s1_idx]
        result[model.s2_idx]  = zdot[model.s2_idx]
        result[model.pq_idx]  = zdot[model.pq_idx]

        return 0

residual_fn = centrifuge_residual()

def solve(model):

    t   = np.empty([model.iterations+1, ], dtype=float)

    WM  = np.empty([model.iterations+1, ], dtype=float)
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
             if model.rb_type == 2:
                 # do regularization for prescribed head on right boundary
                 n_spanning_points = 5

                 z0_view = z0[model.last_idx-n_spanning_points:model.last_idx+1]
                 y_view = model.y[-1-n_spanning_points:]
                 a = ((model.h_last - model.h_init)
                      / (y_view[-1] - y_view[0]) ** 2)

                 z0_view[:] = ((a * np.power(y_view - y_view[0], 2))
                                + model.h_init)
                 z0_view[-1] = model.h_last

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

             WM_total, WM_in_tube = water_mass(u[0, :], mass_in, mass_out,
                                               s1, s2, model)
             WM[0] = WM_total
             if model.calc_gc:
                 GC    = np.empty([model.iterations+1, ], dtype=float)
                 GC[0] = calc_gc(u[0, :], mass_in, mass_out, s1, s2, WM_in_tube,
                                 model)
             if model.calc_rm:
                 RM    = np.empty([model.iterations+1, ], dtype=float)
                 RM[0] = calc_rm(t[0], u[0, :], mass_in, mass_out, s1, s2,
                                 model)
        else:
            z0 = z[i-1, :]

        (flag, t_out, z[i, :]) = simulate_direct(model, residual_fn, z0)

        t[i] = t[i-1] + model.duration

        u[i, :] = h2u(z[i, model.first_idx: model.last_idx+1],
                      model.n, model.m, model.gamma)
        s1 = z[i, model.s1_idx]
        s2 = z[i, model.s2_idx]

        mass_in  = z[i, model.mass_in_idx]
        mass_out = z[i, model.mass_out_idx]
        WM_total, WM_in_tube = water_mass(u[i, :], mass_in, mass_out,
                                          s1, s2, model)
        WM[i] = WM_total
        if model.calc_gc:
            GC[i] = calc_gc(u[i, :], mass_in, mass_out, s1, s2, WM_in_tube,
                            model)
        if model.calc_rm:
            RM[i] = calc_rm(t[i], u[i, :], mass_in, mass_out, s1, s2, model)

        # print('results', GC, z[:, model.mass_out_idx], u,
        #       z[:, :model.last_idx+1])


    if model.draw_graphs:
        from modules.shared.show import draw_graphs

        h = z[:, model.first_idx:model.last_idx+1]

        if not model.calc_gc:
            GC = None

        if not model.calc_rm:
            RM = None

        draw_graphs(t, y=model.y, h=h, u=u, mass_out=z[:, model.mass_out_idx],
                    GC=GC, RM=RM, WM=WM,
                    s1=z[:, model.s1_idx], s2=z[:, model.s2_idx],
                    save_figures=model.save_figures,
                    separate_figures=model.separate_figures)

    if not model.calc_gc:
        GC = np.asarray([], dtype=float)
    if not model.calc_rm:
        RM = np.asarray([], dtype=float)

    return (flag, t, z, GC, RM)
