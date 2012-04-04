import numpy as np

import scikits.odes.sundials.ida as ida
from scikits.odes.sundials.common_defs import ResFunction
from modules.shared.shared_functions import find_omega2g, h2Kh, dudh, h2u

def draw_graphs(fignum, t, x, h, u, mass_out, GC = None, RM = None, WM = None):
    import matplotlib.pyplot as plt

    legend_data = []
    for i in range(len(t)):
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
    plt.plot(t, mass_out, '.')
    plt.xlabel('Time [$s$]')
    plt.ylabel('Outspelled water [$cm^3$]')

    if GC:
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
        q_last  = -Kh_last * np.minimum(0., dhdr_last/ds - omega2g*(r0 + L))
        #print('ql:', q_last)
        #print('q_out', q_last)

        du_dh = dudh(h, n, m, gamma, Ks)
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
                      + np.sum((dy[:-1] + dy[1:])*u[1:-1], 1))

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
        print(model.fl1, model.l0, model.fl2, gc_left)
        gc_right = model.fl1 + model.l0 + model.fl2 - gc_left
        GC       = gc_right
        #print('GC: unsat, sat, gc/wm', gc_unsat, gc_sat, GC[i])
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

    t = np.empty([len(model.duration)+1, ], dtype=float)
    z = np.empty([len(t), model.z_size], dtype=float) # 2 columns: wl_in, wl_out
    u = np.empty([len(t), model.inner_points+2], dtype=float)
    z0  = np.zeros([model.z_size, ], float)
    zp0 = np.zeros([model.z_size, ], float)

    solver = ida.IDA(residual_fn,
                     #compute_initcond='yp0',
                     first_step_size=1e-20,
                     atol=model.atol, rtol=model.rtol,
                     max_step_size=840.,
                     max_steps=8000,
                     #algebraic_vars_idx=[4],
                     linsolver='band', uband=1, lband=1,
                     user_data=model)

    t_end = 0.0

    while model.next_iteration():
        if model.first_iteration_p:
             # initial pressure h
             z0[model.first_idx:model.last_idx+1] = model.h_init
             z0[model.s2_idx] = model.l0

             t[0] = 0.0
             z[0, :] = z0

             solver._init_step(0.0, z0, zp0)

        t_end = t_end + model.duration
        flag, t_out = solver.step(t_end, z[i+1, :])
        t[i+1] = t_out

        u[i+1, :] = h2u(z[i+1, model.first_idx: model.last_idx],
                        model.n, model.m, model.gamma)

        mass_in  = z[:, model.mass_in_idx]
        mass_out = 0.0 # for GC and RM no expelled water is taken into account
        GC[i], _RM, _MW = characteristics(t_out, u, mass_in, mass_out, s1, s2,
                                          model, chtype='gc')


        if t_out < t_end:
            print('Calculation was not finished. Error occured.')
            break

        # falling head
        if model.fh_duration > 0.:
            acc = model.include_acceleration

            model.r0    = model.r0_fall
            model.omega = model.omega_fall
            t_end       = t_end + model.fh_duration

            flag, tfh_out = solver.step(t_end)
            model.include_acceleration = acc

    #print('t,z: ', t, z)

    if model.draw_graphs:
        from shared_functions import y2x

        x = y2x(model.y, z[:, s1_idx], z[:, s2_idx])
        h = z[:, model.first_idx:model.last_idx+1]


        draw_graphs(1, t, x, h, u, z[:, model.mass_out_idx], GC=GC)

    return (flag, t, z)
