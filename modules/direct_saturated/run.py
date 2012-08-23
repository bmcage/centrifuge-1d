import numpy as np

from scikits.odes.sundials.ida import IDA_RhsFunction
from modules.shared.functions import find_omega2g
from modules.shared.solver import simulate_direct

class direct_saturated_rhs(IDA_RhsFunction):
    def evaluate(self, t, x, xdot, result, model):

        omega2g = find_omega2g(t, model)

        L = model.l0
        wl = x[0]
        rS = model.r0 - model.fl1 - wl
        rE = model.r0 + L + model.fl2

        qt = (omega2g/2. / (model.fl1/model.ks1 + L/model.ks
                            + model.fl2/model.ks2)
              * (rE*rE - rS*rS))

        result[model.mass_in_idx]  = xdot[model.mass_in_idx]  + qt
        result[model.mass_out_idx] = xdot[model.mass_out_idx] - qt

        return 0

residual_fn = direct_saturated_rhs()

def solve(model):

    t  = np.empty([model.iterations+1, ], dtype=float)
    z  = np.empty([model.iterations+1, 2], dtype=float) # z[wl_in, wl_out]
    z0 = np.empty([2, ], dtype=float)

    z0[1] = 0.0

    model.init_iteration()

    t[0] = 0.0
    z[0, model.mass_in_idx]  = model.wl0
    z[0, model.mass_out_idx] = 0.0

    while model.next_iteration():
        i = model.iteration
        z0[0] = model.wl0
        #print(i, z0)

        (success_p, t_out, z[i, :]) = simulate_direct(model, residual_fn, z0)
        t[i] = t_out
        #print('zi', z[i, :])
        if not success_p: break

    if model.draw_graphs:
        from modules.shared.show import draw_graphs

        draw_graphs(t, mass_in=z[:, mass_in_idx], mass_out=z[:, mass_out_idx])
    #print('z result:', z[:, :])

    return (success_p, t, z)

def run(model):
    return solve(model)
