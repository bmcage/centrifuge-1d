import numpy as np

from scikits.odes.sundials.ida import IDA_RhsFunction
from modules.shared.functions import right_derivative, y2x
from modules.shared.vangenuchten import h2Kh, dudh, h2u
from modules.direct_draining_saturated.characteristics import \
     water_mass, calc_gc, calc_rm
from modules.shared.solver import simulate_direct
from modules.shared.show import make_dplot, add_dplotline, display_dplots

#TODO: will the new characteristics work also for the previous
#      rb_types?

class centrifuge_residual(IDA_RhsFunction):

    def evaluate(self, t, z, zdot, result, model):
        # F(t,h,dh/dt) = porosity * du/dh * dh/dt
        #                + d[K(h)*(dh/dr - omega^2/g * r)]/dr

        #print('--------------------')
        #print('z', z)
        #print('zdot', zdot)
        #print('--------------------')

        rE = model.re
        L  = model.l0
        fl2 = model.fl2
        rb_type = model.rb_type
        r0 = rE - fl2 - L

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

        omega2g = model.find_omega2g(t)


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


        verbosity = model.verbosity

        if verbosity > 3:
            print('t: %10.6g' % t, 's1: %8.6g' % z[model.s1_idx],
                  's2: %8.6g' % z[model.s2_idx], 'ds2dt: %10.6g' % ds2dt,
                  end='')

        if rb_type == 3:
            q_s2 = -Ks * (dhdy[-1]/ds - omega2g*(r0 + s2))

            rD = rE - model.dip_height
            rI = r0 + s2

            q_sat = (omega2g/2. * (rD*rD - rI*rI)
                     / (model.fl1/model.ks1 + (L-s2)/model.ks
                        + model.fl2/model.ks2))

            if verbosity > 4:
                print('  q_sat', q_sat, 'rD', rD, 'rI', rI,
                      'omega^2/g', omega2g, end='')

            if q_sat < 0:
                if verbosity > 1:
                    print('       Q_SAT: ', q_sat, ' !!!')

                #q_sat = 0.0

            u = h2u(h, n, m, gamma)
            WM_total, WM_in_tube = water_mass(u,
                                              z[model.mass_in_idx],
                                              z[model.mass_out_idx],
                                              z[model.s1_idx], z[model.s2_idx],
                                              model)
            result[last_idx]  = hdot[-1]
            if verbosity > 4:
                print('  WM: ', WM_total, 'WM0', model.wm0, end='')
            result[model.s2_idx] = WM_total - model.wm0
            result[model.mass_out_idx] = zdot[model.mass_out_idx]  - q_sat
        else:
            if rb_type == 0:
                q_out = 0.
                result[last_idx]  = (porosity * du_dh[-1] * hdot[-1]
                                     + 2 / dy[-1] / ds * (q_out - q12[-1]))
            elif rb_type == 1:
                Kh_last =  h2Kh(h[-1], n, m, gamma, Ks)
                q_out  = np.maximum(1e-12,
                                    -Kh_last * (dhdr_last/ds - omega2g*(r0 + L)))
                result[last_idx]  = (porosity * du_dh[-1] * hdot[-1]
                                     + 2 / dy[-1] / ds * (q_out - q12[-1]))
            elif rb_type == 2:
                Kh_last =  h2Kh(h[-1], n, m, gamma, Ks)
                q_out  = np.maximum(1e-12,
                                    -Kh_last * (dhdr_last/ds - omega2g*(r0 + L)))
                result[last_idx]  = hdot[-1]
            else:
                raise NotImplementedError('rb_type has to be 0 - 6')

            result[model.s2_idx]  = zdot[model.s2_idx]
            result[model.mass_out_idx] = zdot[model.mass_out_idx]  - q_out

        result[model.mass_in_idx]  = zdot[model.mass_in_idx]
        result[model.s1_idx]  = zdot[model.s1_idx]
        result[model.pq_idx]  = zdot[model.pq_idx]

        if verbosity > 3: print()

        return 0

residual_fn = centrifuge_residual()

def solve(model):
    # define s1, s2, mass_in, mass_out, u0, wm0, wm_in_tube0
    if model.rb_type == 3:
        if model.s2_0:
            s2 = model.s2_0
        else:
            s2 = model.l0 - model.dip_height
    else:
        s2 = model.l0

    s1 = 0.0
    mass_in = mass_out = 0.0

    wm0 = wm_in_tube0 = 0.0

    u0 = None

    # z0 inicialization
    def initialize_z0(z0, model):
        nonlocal u0, wm0, wm_in_tube0

        z0[model.first_idx:model.last_idx+1] = model.h_init

        if model.rb_type in [2, 3]:
            # do regularization for prescribed head on right boundary
            n_spanning_points = 15

            z0_view = z0[model.last_idx-n_spanning_points:model.last_idx+1]
            y_view = model.y[-1-n_spanning_points:]
            a = ((model.h_last - model.h_init)
                 / (y_view[-1] - y_view[0]) ** 2)

            z0_view[:] = ((a * np.power(y_view - y_view[0], 2))
                      + model.h_init)
            z0_view[-1] = model.h_last

        z0[model.s1_idx] = s1
        z0[model.s2_idx] = s2
        z0[model.mass_in_idx]  = mass_in
        z0[model.mass_out_idx] = mass_out

        # assign value to u0, WM0 and wm0 (wm0 is needed for mass balance)
        u0 = h2u(z0[model.first_idx: model.last_idx+1],
                 model.n, model.m, model.gamma)
        (wm0, wm_in_tube0) = water_mass(u0, mass_in, mass_out, s1, s2, model)
        model.wm0 = wm0

    def initialize_zp0(zp0, z0, model):
        (Ks, n, m, gamma) = (model.ks, model.n, model.m, model.gamma)
        (first_idx, last_idx) = (model.first_idx, model.last_idx)
        h    =  z0[first_idx:last_idx+1]
        h12  = (h[1:] + h[:-1]) / 2
        Kh12 = h2Kh(h12, n, m, gamma, Ks)

        t = 0.0
        omega2g = model.find_omega2g(t)
        (y, dy)  = (model.y, model.dy)
        dhdy12 = (h[1:] - h[:-1]) / dy
        dhdy = np.empty([model.inner_points+2,], dtype=float)
        dhdy[0] = (model.ldc1[0]*h[0] + model.ldc2[0]*h[1] + model.ldc3[0]*h[2])
        dhdy[1:-1] = (model.ldc1[1:-1]*h[:-2] + model.ldc2[1:-1] * h[1:-1]
                     + model.ldc3[1:-1] * h[2:])
        dhdy[-1] = (model.ldc1[-1]*h[-3] + model.ldc2[-1]*h[-2]
                    + model.ldc3[-1]* h[-1])

        s1    = z0[model.s1_idx]
        s2    = z0[model.s2_idx]
        ds    = s2 - s1
        ds1dt = 0.0

        (rE, fl2, L)= (model.re, model.fl2, model.l0)
        r0 = rE - fl2 - L

        q_first = 0.
        q12 = -Kh12 * (dhdy12/ds  - omega2g*(r0 + s1 + ds * model.y12))

        porosity = model.porosity
        du_dh = dudh(h, n, m, gamma)

        rb_type = model.rb_type
        if rb_type == 3:
            (rD, rI) = (rE - model.dip_height, r0 + s2)
            q_sat = (omega2g/2. * (rD*rD - rI*rI)
                     / (model.fl1/model.ks1 + (L-s2)/model.ks
                        + fl2/model.ks2))

            dmodt = q_sat
            zp0_last = 0.0
            ds2dt = +0.1 # some constant for derivation estimate
        else:
            if rb_type == 2:
                zp0_last = 0.0
            else:
                if rb_type == 1:
                    Kh_last =  h2Kh(h[-1], n, m, gamma, Ks)
                    q_out  = np.maximum(1e-12,
                                        -Kh_last*(dhdy_last/ds - omega2g*(r0 + L)))
                else:
                    q_out = 0.0
                zp0_last = \
                  (-2./(porosity * du_dh[-1] * dy[-1]*ds) * (q_out - q12[-1]))
            ds2dt = 0.0

        zp0[last_idx] = zp0_last
        zp0[model.mass_in_idx] = zp0[model.s1_idx] = 0.0
        zp0[model.mass_out_idx] = dmodt
        zp0[model.s2_idx] = ds2dt

        zp0[first_idx] = \
          dhdy[0]/ds*ds1dt - 2./(porosity*du_dh[0]* dy[0]*ds) * (q12[0] - q_first)
        zp0[first_idx+1:last_idx] = \
          (dhdy[1:-1]/ds*((1-y[1:-1])*ds1dt + y[1:-1]*ds2dt)
           - 2./(porosity*du_dh[1:-1]*(dy[:-1] + dy[1:])*ds) * (q12[1:] - q12[:-1]))

    # Initialization
    if model.rb_type == 3:
        atol_backup        = model.atol # backup value
        if type(atol_backup) == list:
            atol = np.asarray(atol_backup, dtype=float)
        else:
            atol = atol_backup * np.ones([model.z_size,], dtype=float)
        atol[model.s2_idx] = model.s2_atol
        model.atol         = atol

        algvars_idx = [model.s2_idx]
    else:
        algvars_idx = None

    if model.estimate_zp0:
        zp0_init = initialize_zp0
    else:
        zp0_init = None

    # Computation
    (flag, t, z) = simulate_direct(initialize_z0, model, residual_fn,
                                   root_fn = None, nr_rootfns=None,
                                   initialize_zp0=zp0_init,
                                   algvars_idx=algvars_idx)

    # Restore modified values
    model.atol = atol_backup

    # Results
    k  = np.alen(t)
    s1 = z[:, model.s1_idx]
    s2 = z[:, model.s2_idx]
    mass_in  = z[:, model.mass_in_idx]
    mass_out = z[:, model.mass_out_idx]

    no_measurements = np.empty([0,], dtype=float)

    if model.calc_wm:
        h  = z[:, model.first_idx: model.last_idx+1]
        u  = np.empty([k, model.inner_points+2], dtype=float)

        WM = np.empty(t.shape, dtype=float)
        WM_in_tube = np.empty(t.shape, dtype=float)

        u[0, :] = u0
        WM[0]   = wm0
        WM_in_tube[0] = wm_in_tube0
        for i in range(1, k):
            u[i, :] = h2u(h[i, :], model.n, model.m, model.gamma)

            (wm_total, wm_in_tube) = \
              water_mass(u[i, :], mass_in[i], mass_out[i], s1[i], s2[i], model)
            WM[i]         = wm_total
            WM_in_tube[i] = wm_in_tube
    else:
        u          = no_measurements
        WM         = no_measurements
        WM_in_tube = no_measurements

    if model.calc_gc:
        GC = np.empty(t.shape, dtype=float)

        for i in range(k):
            GC[i] = calc_gc(u[i, :], mass_in[i], s1[i], s2[i], WM_in_tube[i],
                            model)
    else:
        GC = no_measurements

    if model.calc_rm:
        RM = np.empty(t.shape, dtype=float)

        for i in range(k):
            RM[i] = calc_rm(t[i], u[i, :], mass_in[i], s1[i], s2[i], model)
    else:
        RM = no_measurements

    return (flag, t, z, GC, RM, u, WM, WM_in_tube)

def get_refencing_models(model):

    def models_generator(ref_params):
        if type(ref_params) == dict: # single reference
            ref_params = [ref_params]

        for ref in ref_params:
            backup_params = {}
            for (key, value) in ref.items(): # backup
                if key in model._iterable_parameters:
                    print('Referencing model cannot have different iterable '
                          'parameters than original model:', key)
                    exit(1)
                backup_params[key] = getattr(model, key)

            model.set_parameters(ref)

            yield model

            model.set_parameters(backup_params) # restore

    if not model.params_ref:
        return None
    else:
        return models_generator(model.params_ref)


def multiple_solves(c_model, referencing_models=[]):
    def iterate_models():
        yield c_model

        if referencing_models:
            for model in referencing_models:
                yield model
        else:
            yield None

    collected_computations = []

    for model in iterate_models():

        if not  model: break

        (flag, t, z, GC, RM, u, WM, WM_in_tube) = solve(model)

        if not flag:
            print('For given model the solver did not find results. Skipping.')
            continue
        else:
            any_data = True

        s1 = z[:, model.s1_idx]
        s2 = z[:, model.s2_idx]
        x = y2x(model.y, s1, s2).transpose()
        h = z[:, model.first_idx:model.last_idx+1].transpose()
        u = u.transpose()
        MO = z[:, model.mass_out_idx]
        MI = z[:, model.mass_in_idx]
        collected_computations.append(((t, h, u, GC, RM, WM, MI, MO, s1, s2, x)))

    data_annotation = ('t', 'h', 'u', 'GC', 'RM', 'WM', 'MI', 'MO',
                       's1', 's2', 'x')

    return (collected_computations, data_annotation)

def display_graphs(model, computations, annotation, options):
    from collections import defaultdict

    if not computations: return

    dplots_names = ['h', 'u', 'GC', 'RM', 'WM', 'MI', 'MO', 's1', 's2']
    dplots_bucket = {name: make_dplot(name, legend_loc=1, legend_title=None)
                     for name in dplots_names}
    for name in ['h', 'u']:
        if (not model.separate_figures) and (name == 'h'):
            dplots_bucket[name]['show_legend'] = False
        dplots_bucket[name]['legend_title'] = 'Time [min]'
        dplots_bucket[name]['legend_bbox'] = (1.02, 1.)
        dplots_bucket[name]['legend_loc'] = 2

    dplots_bucket['MO']['legend_loc'] = 4

    line_label = 'computed'

    for (idx, computed_data) in enumerate(computations):
        measurement = dict(zip(annotation, computed_data))

        t = [ti/60. for ti in measurement['t']] # sec -> min
        t_legend = ['% 7d' % ti for ti in t]

        if idx == 0:
            # dirty hack, for now we take as 't' for measurements for models[0]
            # of the first model
            t_meas = t

            # ok, h and u we display only from the first model
            for m in ['h', 'u']:
                add_dplotline(dplots_bucket[m], measurement['x'],
                              measurement[m],
                              label=t_legend, line_opts='-')

        for m in dplots_names[2:]:
            add_dplotline(dplots_bucket[m], t, measurement[m],
                          label=line_label, line_opts='.')

        if idx == 0:
            # and the subsequent lines will be labeled as Reference (default)
            line_label = None

    # Now include also measurements
    t1_meas = t_meas[1:]
    add_dplotline(dplots_bucket['GC'], t1_meas, model.gc1,
                  label='measured', line_opts='x')
    add_dplotline(dplots_bucket['RM'], t1_meas, model.rm1,
                  label='measured', line_opts='x')
    add_dplotline(dplots_bucket['MI'], t1_meas, model.wl1,
                  label='measured', line_opts='x')
    if model.wl_out:
        add_dplotline(dplots_bucket['MO'], t1_meas, np.cumsum(model.wl_out),
                      label='measured', line_opts='x')

    # put it together and display
    dplots = list(dplots_bucket.values())
    display_dplots(dplots, save_figures=options['save_figures'],
                   separate_figures=options['separate_figures'],
                   save_as_text=options['save_as_text'],
                   show_figures=options['show_figures'],
                   experiment_info=options['experiment_info'])


def run(model):
    referencing_models = get_refencing_models(model)
    (results, annotation) = multiple_solves(model, referencing_models)
    display_options = {'save_figures': model.save_figures,
                       'separate_figures': model.separate_figures,
                       'save_as_text': model.save_as_text,
                       'show_figures': model.show_figures,
                       'experiment_info': model.experiment_information}
    display_graphs(model, results, annotation, display_options)
