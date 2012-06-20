import numpy as np

from scikits.odes.sundials.common_defs import IDA_RhsFunction
from modules.shared.shared_functions import find_omega2g, right_derivative
from modules.shared.vangenuchten import h2Kh, dudh, h2u
from modules.direct_draining_saturated.characteristics import \
     water_mass, calc_gc, calc_rm
from modules.shared.solver import simulate_direct

#TODO: will the new characteristics work also for the previous
#      rb_types?

class centrifuge_residual(IDA_RhsFunction):

    def evaluate(self, t, z, zdot, result, model):
        # F(t,h,dh/dt) = porosity * du/dh * dh/dt
        #                + d[K(h)*(dh/dr - omega^2/g * r)]/dr

        #print('--------------------')
        print('t', t)
        #print('z', z)
        #print('zdot', zdot)
        #print('--------------------')

        r0 = model.r0
        L  = model.l0
        rb_type = model.rb_type

        s2 = z[model.s2_idx]
        s1 = z[model.s1_idx]
        ds = s2 - s1

        if (rb_type > 3) and (s2 > L - model.dip_height):
            print('refine, s2EQ=', L - model.dip_height, 's2=', s2)
            return 1

        first_idx = model.first_idx
        last_idx  = model.last_idx

        h    =  z[first_idx:last_idx+1]

        if np.any(h > 0): # positive pressure - we want to refine the step
            return 1
        #print('h: ', h)

        dhdy = np.empty([model.inner_points+2,], dtype=float)

        hdot =  zdot[first_idx:last_idx+1]
        h12  = (h[1:] + h[:-1]) / 2

        omega2g = find_omega2g(t, model)


        Ks = model.ks
        n  = model.n
        m  = model.m
        gamma = model.gamma

        porosity = model.porosity

        Kh12 = h2Kh(h12, n, m, gamma, Ks)

        #Kh12 = cProfile.run('h2Kh(h12, n, m, gamma, Ks)', 'h2kh1')
        #Kh_last =  cProfile.run('h2Kh(h[-1], n, m, gamma, Ks)', 'h2kh2')

        y  = model.y
        dy = model.dy

        dhdy12 = (h[1:] - h[:-1]) / dy
        dhdy[0] = (model.ldc1[0] * h[0]
                     + model.ldc2[0] * h[1]
                     + model.ldc3[0] * h[2])
        dhdy[1:-1] = (model.ldc1[1:-1] * h[:-2]
                     + model.ldc2[1:-1] * h[1:-1]
                     + model.ldc3[1:-1] * h[2:])
        dhdy[-1] = (model.ldc1[-1] * h[-3]
                     + model.ldc2[-1] * h[-2]
                     + model.ldc3[-1] * h[-1])

        #dhdy[:-1] = 0.
        #print('dhdy', dhdy, 'h[1:3]', h[0:3])
        #print('ldc', model.ldc1[:10], model.ldc2[:10], model.ldc3[:10])
        q_first = 0.
        q12 = -Kh12 * (dhdy12/ds  - omega2g*(r0 + s1 + ds * model.y12))

        #print('ql:', q_last)
        #print('q_out', q_last)
        ds1dt = zdot[model.s1_idx]
        ds2dt = zdot[model.s2_idx]

        du_dh = dudh(h, n, m, gamma)
        result[first_idx] = \
          (porosity * du_dh[0] * (hdot[0] - dhdy[0]/ds*ds1dt)
           + 2 / dy[0] / ds * (q12[0] - q_first))

        result[first_idx+1:last_idx] = \
          (porosity*du_dh[1:-1]*(hdot[1:-1]
                                 - dhdy[1:-1]/ds*((1-y[1:-1])*ds1dt
                                                  + y[1:-1]*ds2dt))
            + 2 / (dy[:-1] + dy[1:]) / ds * (q12[1:] - q12[:-1]))

        #TODO: in q_s2 h<0 => K(h) != Ks
        q_s2 = -Ks * (dhdy[-1]/ds - omega2g*(r0 + s2))

        rD = model.r0 + L - model.dip_height
        rI = model.r0 + s2

        q_sat = (omega2g/2. * (rD*rD - rI*rI)
                           / (model.fl1/model.ks1 + (L-s2)/model.ks
                              + model.fl2/model.ks2))
        q_sat = 0.0
        if q_sat < 0:
            print(10*'-' + '\nQ_SAT =', q_sat, ' !!!\n' + 10*'-')
            print(omega2g, rD, rI, L, s2, model.ks, model.ks2, model.ks1)

        print('h3:', h[-3:])
        print('s1', z[model.s1_idx], 's2', z[model.s2_idx])

        if rb_type == 3:
            if dhdy[-1] == 0.0:
                dhdy[-1] = -1e-12
            q_last = q_sat - q_s2
            #q_last = q_sat - q_s2
            dhdt_last = - q_last/ porosity / du_dh[-1]
            dqdy_last = right_derivative([(dy[-2]+dy[-1])/2, dy[-1]/2],
                                         [q12[-2], q12[-1], q_s2])

            #print('s2', z[model.s2_idx], '\ns2dt', -dhdt_last/dhdy[-1])
            print('\nq_sat', q_sat, 'q_s2', q_s2, 'q_last', q_last,
                  'expelled', z[model.mass_out_idx])

            result[last_idx]  = hdot[-1]
            result[model.s2_idx] = ds2dt + dqdy_last/dhdy[-1]/porosity/du_dh[-1]
        elif rb_type == 4:
            q_last = q_s2 - q_sat

            result[last_idx]  = hdot[-1]
            result[model.s2_idx]  = \
              zdot[model.s2_idx] + q_last/porosity
            print('\nq_sat', q_sat, 'q_s2', q_s2, 'q_last', q_last,
                  'expelled', z[model.mass_out_idx])
        elif rb_type == 5:
            rE = model.r0 + L + model.fl2
            rI = model.r0 + s2

            q_out = (omega2g/2. / (model.fl1/model.ks1 + L/model.ks
                                   + model.fl2/model.ks2)
                     * (rE*rE - rI*rI))

            wdot = ds/2  * (dy[0]* hdot[0] + dy[-1]*hdot[-1]
                            + np.sum((dy[:-1] + dy[1:])*hdot[1:-1]))
        else:
            raise NotImplementedError('rb_type has to be 3 - 5')

        #print('ds2dt', result[model.s2_idx], dhdotdr_last/dhdr[-1],
        #      dhdotdr_last, dhdr[-1])

        result[model.mass_in_idx]  = zdot[model.mass_in_idx]
        result[model.mass_out_idx] = zdot[model.mass_out_idx]  - q_sat
        result[model.s1_idx]  = zdot[model.s1_idx]
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

    # initialization:
    # set values for t[0], z[0], u[0], GC[0]
    z0 = z[0, :]

    if model.rb_type > 3:
        atol_orig = model.atol
        atol = atol_orig * np.ones([model.z_size,], dtype=float)
        atol[model.s2_idx] = 1e-4

        rtol_orig = model.rtol
        rtol = rtol_orig * np.ones([model.z_size,], dtype=float)
        rtol[model.s2_idx] = 1e-4

    z0[model.first_idx:model.last_idx+1] = model.h_init
    if model.rb_type in [2, 3, 4]:
        # do regularization for prescribed head on right boundary
        n_spanning_points = 15

        z0_view = z0[model.last_idx-n_spanning_points:model.last_idx+1]
        y_view = model.y[-1-n_spanning_points:]
        a = ((model.h_last - model.h_init)
             / (y_view[-1] - y_view[0]) ** 2)

        z0_view[:] = ((a * np.power(y_view - y_view[0], 2))
                      + model.h_init)
        z0_view[-1] = model.h_last

    if model.rb_type < 3:
        s2 = model.l0
    else:
        #s2 = model.l0 - model.dip_height
        s2 = 1.00

    s1 = 0.0
    mass_in = mass_out = 0.0

    z0[model.s1_idx] = s1
    z0[model.s2_idx] = s2
    z0[model.mass_in_idx]  = mass_in
    z0[model.mass_out_idx] = mass_out

    t[0] = 0.0
    u[0, :] = h2u(z0[model.first_idx: model.last_idx+1],
                  model.n, model.m, model.gamma)

    WM_total, WM_in_tube = water_mass(u[0, :], mass_in, mass_out,
                                      s1, s2, model)
    WM[0] = WM_total
    if model.calc_gc:
        GC    = np.empty([model.iterations+1, ], dtype=float)
        GC[0] = calc_gc(u[0, :], mass_in, s1, s2, WM_in_tube,
                                 model)
    if model.calc_rm:
        RM    = np.empty([model.iterations+1, ], dtype=float)
        RM[0] = calc_rm(t[0], u[0, :], mass_in, mass_out, s1, s2,
                        model)

    while model.next_iteration():
        i = model.iteration

        z0 = z[i-1, :]

        (flag, t_out, z[i, :]) = simulate_direct(model, residual_fn, z0)

        t[i] = t[i-1] + t_out

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
            GC[i] = calc_gc(u[i, :], mass_in, s1, s2, WM_in_tube,
                            model)

        if model.calc_rm:
            RM[i] = calc_rm(t[i], u[i, :], mass_in, mass_out, s1, s2, model)

        # print('results', GC, z[:, model.mass_out_idx], u,
        #       z[:, :model.last_idx+1])

    if model.rb_type > 3:
        model.atol = atol_orig
        model.rtol = rtol_orig

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
                    separate_figures=model.separate_figures,
                    save_as_text=model.save_as_text,
                    model=model)

    if not model.calc_gc:
        GC = np.asarray([], dtype=float)
    if not model.calc_rm:
        RM = np.asarray([], dtype=float)

    return (flag, t, z, GC, RM)
