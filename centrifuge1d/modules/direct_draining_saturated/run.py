from __future__ import division, print_function

import numpy as np

from scikits.odes.sundials.ida import IDA_RhsFunction
from ..shared.functions import y2x
from ..shared.solver import simulate_direct
from ..shared.show import show_results

#TODO: will the new characteristics work also for the previous
#      rb_types?

class centrifuge_residual(IDA_RhsFunction):

    def evaluate(self, t, z, zdot, result, model):
        # F(t,h,dh/dt) = porosity * du/dh * dh/dt
        #                + d[K(h)*(dh/dr - omega^2/g * r)]/dr

        (rE, L, fl2) = (model.re, model.l0, model.fl2)
        r0 = rE - fl2 - L

        (s1, s2) = (z[model.s1_idx], z[model.s2_idx])
        ds = s2 - s1

        rb_type = model.rb_type

        if (rb_type > 3) and (s2 > L - model.dip_height):
            ## TODO: ERROR, no >3 and should be r0+s2 ????
            print('refine, s2EQ=', L - model.dip_height, 's2=', s2)
            return 1

        (first_idx, last_idx) = (model.first_idx, model.last_idx)
        SC = model.SC

        h    =  z[first_idx:last_idx+1]

        if np.any(h > 0): # positive pressure - we want to refine the step
            return 1

        hdot =  zdot[first_idx:last_idx+1]

        Ks = model.ks
        h12  = (h[1:] + h[:-1]) / 2
        Kh12 = SC.h2Kh(h12, Ks)

        (y, dy)  = (model.y, model.dy)

        (ldc1, ldc2, ldc3) = (model.ldc1, model.ldc2, model.ldc3)
        dhdy12 = (h[1:] - h[:-1]) / dy
        dhdy = np.empty([model.inner_points+2,], dtype=float)
        dhdy[0]    = ldc1[0] * h[0] + ldc2[0] * h[1] + ldc3[0]*h[2]
        dhdy[1:-1] = ldc1[1:-1]*h[:-2] + ldc2[1:-1]*h[1:-1] + ldc3[1:-1]*h[2:]
        dhdy[-1]   = ldc1[-1] * h[-3] + ldc2[-1] * h[-2] + ldc3[-1] * h[-1]

        omega2g = model.find_omega2g(t)

        q_first = 0.
        q12 = -Kh12 * (dhdy12/ds  - omega2g*(r0 + s1 + ds * model.y12))

        (ds1dt, ds2dt) = (zdot[model.s1_idx], zdot[model.s2_idx])

        porosity = model.porosity
        du_dh = SC.dudh(h)

        result[first_idx] = \
          (porosity * du_dh[0] * (hdot[0] + dhdy[0]/ds*ds1dt)
           + 2 / dy[0] / ds * (q12[0] - q_first))

        result[first_idx+1:last_idx] = \
          (porosity*du_dh[1:-1]*(hdot[1:-1]
                                 + dhdy[1:-1]/ds*((1-y[1:-1])*ds1dt
                                                  + y[1:-1]*ds2dt))
            + 2 / (dy[:-1] + dy[1:]) / ds * (q12[1:] - q12[:-1]))

        if rb_type == 3:
            rD = rE - model.dip_height
            rI = r0 + s2

            #no fl1 contribution, as no saturated flow there
            q_sat = (omega2g/2. * (rD*rD - rI*rI)
                     / ( (L-s2)/model.ks
                        + model.fl2/model.ks2))

            if q_sat < 0:
                pass
                #q_sat = 0.0

            u = SC.h2u(h)
            (WM_total, WM_in_tube) = \
              model.measurements.store_calc_wm(u, model.dy, s1, s2,
                                               z[model.mass_in_idx], L-s2,
                                               z[model.mass_out_idx],
                                               model.porosity, model.fl2,
                                               model.fp2, store=False)
            result[last_idx]  = hdot[-1]
            result[model.s2_idx] = WM_total - model.wm0
            result[model.mass_out_idx] = zdot[model.mass_out_idx]  - q_sat
        else:
            if rb_type == 0:
                q_out = 0.
                result[last_idx]  = (porosity * du_dh[-1] * hdot[-1]
                                     + 2 / dy[-1] / ds * (q_out - q12[-1]))
            elif rb_type == 1:
                Kh_last =  SC.h2Kh(h[-1], Ks)
                q_out  = np.maximum(1e-12,
                                    -Kh_last * (dhdy[-1]/ds - omega2g*(r0 + L)))
                result[last_idx]  = (porosity * du_dh[-1] * hdot[-1]
                                     + 2 / dy[-1] / ds * (q_out - q12[-1]))
            elif rb_type == 2:
                Kh_last =  SC.h2Kh(h[-1], Ks)
                q_out  = np.maximum(1e-12,
                                    -Kh_last * (dhdy[-1]/ds - omega2g*(r0 + L)))
                result[last_idx]  = hdot[-1]
            else:
                raise NotImplementedError('rb_type has to be 0 - 3')

            result[model.s2_idx]  = zdot[model.s2_idx]
            result[model.mass_out_idx] = zdot[model.mass_out_idx]  - q_out

        result[model.mass_in_idx] = zdot[model.mass_in_idx]
        result[model.s1_idx]      = zdot[model.s1_idx]
        result[model.pq_idx]      = zdot[model.pq_idx]

        verbosity = model.verbosity
        if verbosity > 2:
            if (rb_type == 3) and (q_sat < 0.0):
                print('       Q_SAT: ', q_sat, ' !!!')
            if verbosity > 3:
                print('t: %10.6g' % t, 's1: %8.6g' % z[model.s1_idx],
                      's2: %8.6g' % z[model.s2_idx], 'ds2dt: %10.6g' % ds2dt)

                if (verbosity > 4) and (rb_type == 3):
                    print('  q_sat', q_sat, 'rD', rD, 'rI', rI,
                          'omega^2/g', omega2g,
                          '  WM: ', WM_total, 'WM0', model.wm0)

        return 0

RESIDUAL_FN = centrifuge_residual()

def on_measurement(t, z, model, measurements):

    s1 = measurements.store_calc_measurement('s1', z[model.s1_idx])
    s2 = measurements.store_calc_measurement('s2', z[model.s2_idx])
    MI = measurements.store_calc_measurement('MI', z[model.mass_in_idx])
    MO = measurements.store_calc_measurement('MO', z[model.mass_out_idx])
    L  = measurements.store_calc_measurement('L', model.l0) #height sample is fixed
    x  = y2x(model.y, s1, s2)

    u = measurements.store_calc_u(x, z[model.first_idx: model.last_idx+1],
                                  model.SC)

    (WM, WM_in_tube) = \
      measurements.store_calc_wm(u, model.dy, s1, s2, MI,
                                 model.l0-s2, MO, model.porosity,
                                 model.fl2, model.fp2)

    if measurements.calc_measurement_p('GC'):
        measurements.store_calc_gc(u, model.y, model.dy, s1, s2,
                                   MI, s2, model.l0, model.porosity,
                                   model.fl2, model.fp2, model.l0,
                                   WM_in_tube, model.density,
                                   from_end=model.l0 + model.fl2)

    if measurements.calc_measurement_p('RM'):
        measurements.store_calc_rm(t, u, MI, s1, s2, model)

    omega2g = model.find_omega2g(t)

    if measurements.calc_measurement_p('gF_MT'):
        l0 = model.l0
        r0 = model.re - model.fl2 - l0

        measurements.store_calc_gf_mt(omega2g, u, model.y, model.dy, r0,
                                      s1, s2, MI, s2, l0, model.porosity,
                                      model.fl2, model.fp2, l0, l0,
                                      model.density,
                                      model.tube_crosssectional_area)

    if measurements.calc_measurement_p('gF_MO'):
        measurements.store_calc_gf_mo(omega2g, MO,
                                      model.mo_gc_calibration_curve,
                                      model.density,
                                      model.tube_crosssectional_area)

def initialize_z0(z0, model):
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

    if model.rb_type == 3:
        if model.s2_0:
            s2 = model.s2_0
        else:
            s2 = model.l0 - model.dip_height
    else:
        s2 = model.l0

    s1 = 0.0
    mass_in = mass_out = 0.0

    z0[model.s1_idx] = s1
    z0[model.s2_idx] = s2
    z0[model.mass_in_idx]  = mass_in
    z0[model.mass_out_idx] = mass_out

    # assign value to u0, WM0 and wm0 (wm0 is needed for mass balance)
    u0 = model.SC.h2u(z0[model.first_idx: model.last_idx+1])

    (wm0, wm_in_tube0) = \
      model.measurements.store_calc_wm(u0, model.dy, s1, s2, mass_in,
                                       model.l0-s2, mass_out, model.porosity,
                                       model.fl2, model.fp2, store=False)

    model.wm0 = wm0

def initialize_zp0(zp0, z0, model):
    Ks = model.ks
    (first_idx, last_idx) = (model.first_idx, model.last_idx)
    h    =  z0[first_idx:last_idx+1]
    h12  = (h[1:] + h[:-1]) / 2
    Kh12 = model.SC.h2Kh(h12, Ks)

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
    du_dh = model.SC.dudh(h)

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
            Kh_last =  model.SC.h2Kh(h[-1], Ks)
            q_out  = np.maximum(1e-12,
                                -Kh_last * (dhdy[-1]/ds - omega2g*(r0 + L)))
        else:
            if rb_type == 1:
                Kh_last =  model.SC.h2Kh(h[-1], Ks)
                q_out  = np.maximum(1e-12,
                                    -Kh_last*(dhdy[-1]/ds - omega2g*(r0 + L)))
            else:
                q_out = 0.0

            zp0_last = \
              (-2./(porosity * du_dh[-1] * dy[-1]*ds) * (q_out - q12[-1]))

        dmodt = q_out
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

def solve(model, measurements):
    measurements.reset_calc_measurements()

    if model.dynamic_h_init:
        model.h_init = model.SC.get_dyn_h_init(model.c_gammah, model.h_init_max)
        if model.verbosity > 2:
            print('\nh_init: ', model.h_init,
            'u_init', model.SC.h2u(model.h_init), '\n')

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
    (flag, t, z, i) = simulate_direct(initialize_z0, model, model.measurements,
                                      RESIDUAL_FN,
                                      root_fn = None, nr_rootfns=None,
                                      initialize_zp0=zp0_init,
                                      algvars_idx=algvars_idx,
                                      on_measurement=on_measurement)

    # Restore modified values
    if model.rb_type == 3:
        model.atol = atol_backup

    return flag

solve_direct = solve  # to remove warning about 'solve_direct' not specified

import cProfile, pstats
PROFILE = False

def run(model):
    if PROFILE:
        pr = cProfile.Profile()
        pr.enable()

    show_results(model.experiment_info, model=model)

    if PROFILE:
        pr.disable()
        pr.print_stats('cumtime')

def dry_run(model):
    show_results(model.experiment_info, model=model)
