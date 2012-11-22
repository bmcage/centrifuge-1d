from __future__ import division

import numpy as np

from scikits.odes.sundials.ida import IDA_RhsFunction
from modules.shared.solver import simulate_direct
from modules.shared.characteristics import calc_f_mo

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

def solve(model):

    def initialize_z0(z0, model):
        z0[model.mass_in_idx]  = model.wl0
        z0[model.mass_out_idx] = 0.0

    def update_init(i, z0, model):
        z0[model.mass_in_idx] = model.wl0

    def add_extra_measurement(t, z, model, measurements):
         measurements['omega2g'].append(model.find_omega2g(t))

    measurements = {'omega2g': []}

    (flag, t, z) = simulate_direct(initialize_z0, model, residual_fn,
                                   update_initial_condition=update_init,
                                   measurements=measurements,
                                   take_measurement=add_extra_measurement)

    MI = z[:, model.mass_in_idx].transpose()
    MO = z[:, model.mass_out_idx].transpose()

    measurements['MI'] = MI
    measurements['MO'] = MO

    if model.calc_f_mo:
        F_MO = np.empty(MO.shape, dtype=float)

        omega2gs = measurements['omega2g']

        for (i, omega2g) in enumerate(omega2gs):
            F_MO[i] = calc_f_mo(omega2g, MO[i], model.mo_gc_calibration_curve)

        measurements['F_MO'] = F_MO

    return (flag, t, z, measurements)

def extract_data(model):
    (flag, t, z, measurements) = solve(model)

    if not flag:
        print('For given model the solver did not find results. Skipping.')

    extracted_data = {name: (t, value) for (name, value) in measurements}

    return (flag, extracted_data)

def run(model):
    from modules.shared.functions import show_results
    show_results(model.experiment_info, model=model)
