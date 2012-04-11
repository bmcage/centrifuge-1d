import numpy as np

from scikits.odes.sundials.common_defs import ResFunction
from modules.shared.shared_functions import find_omega2g
from modules.shared.solver import simulate_direct

def draw_graphs(fignum, t, wl_in, wl_out):

    import matplotlib.pyplot as plt

    rows = 1
    figheight = 4

    plt.figure(fignum, figsize=(12, figheight))

    plt.subplot(rows,2,1)
    plt.plot(t, wl_in, 'b')
    plt.xlabel('Time')
    plt.ylabel('Water in inflow chamber [cm]')

    plt.subplot(rows,2,2)
    plt.plot(t, wl_out, 'k')
    plt.xlabel('Time')
    plt.ylabel('Water in outflow chamber [cm]')

    plt.show()

class direct_saturated_rhs(ResFunction):
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

    t = np.empty([model.iterations+1, ], dtype=float)
    z = np.empty([model.iterations+1, 2], dtype=float) # z[wl_in, wl_out]

    model.init_iteration()
    t[0] = 0.0
    z[0, :] = np.asarray([model.wl0, 0.0], dtype = float)

    while True:
        i = model.iteration
        z0 = z[i-1, :]

        (success_p, t_out, z[i, :]) = simulate_direct(model, residual_fn, z0)

        t[i] = t_out

        if not (success_p and model.next_iteration()): break

    if model.draw_graphs:
        draw_graphs(1, t, z[:, mass_in_idx], z[:, mass_out_idx])

    return (success_p, t, z)
