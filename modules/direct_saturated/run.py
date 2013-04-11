from __future__ import division

import numpy as np

from scikits.odes.sundials.ida import IDA_RhsFunction
from modules.shared.solver import simulate_direct

class direct_saturated_rhs(IDA_RhsFunction):
    def evaluate(self, t, x, xdot, result, model):

        omega2g = model.find_omega2g(t)

        wl  = x[0]
        rE  = model.re
        L   = model.l0
        fl2 = model.fl2
        r0  = rE - fl2 - L
        rS  = r0 - model.fl1 - wl

        qt = (omega2g/2. / (model.fl1/model.ks1 + L/model.ks
                            + fl2/model.ks2)
              * (rE*rE - rS*rS))

        result[model.mass_in_idx]  = xdot[model.mass_in_idx]  + qt
        result[model.mass_out_idx] = xdot[model.mass_out_idx] - qt

        return 0

residual_fn = direct_saturated_rhs()

def on_measurement(t, z, model, measurements):
    MI = measurements.store_calc_measurement('MI', z[model.mass_in_idx])
    MO = measurements.store_calc_measurement('MO', z[model.mass_out_idx])

    if measurements.calc_measurement_p('gF_MO'):
        omega2g = model.find_omega2g(t)

        measurements.store_calc_gf_mo(omega2g, MO,
                                      model.mo_gc_calibration_curve,
                                      model.density,
                                      model.tube_crosssectional_area)

    if measurements.calc_measurement_p('gF_MT'):
        print("The 'gf_mt' option is not supported yet for direct saturated "
              "flow. Skipping.")
        #exit(1)

def initialize_z0(z0, model):
    z0[model.mass_in_idx]  = model.wl0
    z0[model.mass_out_idx] = 0.0

    print('z0', z0)

def update_init(i, z0, model):
    z0[model.mass_in_idx] = model.wl0

def solve(model, measurements):
    (flag, t, z, i) = \
      simulate_direct(initialize_z0, model, measurements, residual_fn,
                      update_initial_condition=update_init,
                      on_measurement=on_measurement)

    return flag

def run(model):
    from modules.shared.show import show_results
    show_results(model.experiment_info, model=model)
