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

def solve(model):

    def initialize_z0(z0, model):
        z0[model.mass_in_idx]  = model.wl0
        z0[model.mass_out_idx] = 0.0

    def update_init(i, z0, model):
        z0[model.mass_in_idx] = model.wl0

    (flag, t, z) = simulate_direct(initialize_z0, model, residual_fn,
                                   update_initial_condition=update_init)

    return (flag, t, z)

def extract_data(model):
    (flag, t, z) = solve(model)

    if not flag:
        print('For given model the solver did not find results. Skipping.')

    extracted_data = {'MI': (t, z[:, model.mass_in_idx]),
                      'MO': (t, z[:, model.mass_out_idx])}
    return (flag, extracted_data)

def run(model):
    from modules.shared.functions import show_results
    show_results(model)
