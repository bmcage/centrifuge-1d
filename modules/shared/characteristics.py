import numpy as np

def water_mass(u, mass_in, mass_out, s1, s2, model):
    """
      Determine the amount of water contained in the experiment.
      The amount is the sum of water in inflow chamber, outflow chamber,
      saturated zone in <r0, s1> and unsaturates zone in <s1, s2>.
    """
    dy = model.dy
    ds = s2 - s1

    # Water mass
    wm_sat   = s1
    wm_unsat = ds/2  * (dy[0]* u[0] + dy[-1]*u[-1]
                        + np.sum((dy[:-1] + dy[1:])*u[1:-1]))

    WM_in_tube = model.density * (mass_in + model.porosity*(wm_sat + wm_unsat))
    WM_total   = WM_in_tube + model.density * mass_out

    return WM_total, WM_in_tube

def calc_gc(u, mass_in, mass_out, s1, s2, WM_in_tube, model):
    """
      Determine the gravitational center in the sample.
      GC is found from the start of the soil sample (not from centr.axis)
    """
    y  = model.y
    dy = model.dy
    ds = s2 - s1
    L  = model.l0
    l0_out = L + model.wt_out
    l_out  = L + model.wt_out - mass_out

    # sample = filter1 + soil + filter2
    r0 = model.fl1

    gc_unsat = (model.porosity * 1/2 * model.density * ds
                * ((r0 + s1)*dy[0]*u[0]
                   + (r0 + s2)*dy[-1]*u[-1]
                   + np.sum((dy[:-1]+dy[1:])
                            *(r0 + s1 + ds*y[1:-1])*u[1:-1])))
    gc_sat   = (1/2 * model.density
                * (model.porosity * (np.power(r0 + s1, 2) - np.power(r0, 2))
                   + (np.power(r0, 2) - np.power(r0 - mass_in, 2))
                   + (np.power(r0 + l0_out, 2) - np.power(r0 + l_out, 2))))
    gc_from_left  = (gc_unsat + gc_sat) / WM_in_tube
    gc_from_right = model.fl1 + model.l0 + model.fl2 - gc_from_left

    GC            = gc_from_right

    return GC

def calc_rm(t, u, mass_in, mass_out, s1, s2, model):

    porosity = model.porosity
    y  = model.y
    dy = model.dy
    L  = model.l0
    l0_out = L + model.wt_out
    l_out  = L + model.wt_out - mass_out

    ds = s2 - s1

    P = np.pi * model.d / 4
    omega2g = find_omega2g(t, model.omega, model)

    # Rotational momentum
    r0 = model.r0
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

    return RM
