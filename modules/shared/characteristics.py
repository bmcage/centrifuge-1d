import numpy as np

def water_mass(u, dy, s1, s2, mass_in, saturated_length, free_fluid_length,
               porosity, fl2, fp2):
    """
      Determine the amount of water contained in the experiment.
      The amount is the sum of water in inflow chamber (mass_in), outflow
      chamber + any other source like water in basin (free_fluid_length)
      and saturated zone (saturated_length) just as unsaturated zone
      (which is between s1 and s2). The soil is characterized by it's
      porosity (porosity).
      Also water in filter of length fl2 with porosity fp2 is considered,
      but only for fully saturated filter.
    """
    # Water mass
    ds = s2 - s1
    unsat = ds/2  * (dy[0]* u[0] + dy[-1]*u[-1]
                     + np.sum((dy[:-1] + dy[1:])*u[1:-1]))

    WM_in_tube = (mass_in + porosity*(saturated_length + unsat))
    WM_in_tube += fp2 * fl2

    WM_total   = WM_in_tube + free_fluid_length

    return WM_total, WM_in_tube

def calc_gc(u, mass_in, s1, s2, WM_in_tube, model):
    """
      Determine the gravitational center of water in the sample. The water
      on thi inflow is taken into account (but not water on the outflow).
      Computed GC is measured from the end of the soil sample
      (in the direction from centrifuge axis)
    """

    if model.fl1 > 0.0:
        raise ValueError('GC not implemented for starting filter.')

    L   = model.l0
    fl2 = model.fl2

    # We assume that ending filter remains fully saturated, so the case for
    # unsaturated (i.e. sample at the end has u<0.999) is not implemented
    if (fl2 > 0.0) and (s2 > (L - 1e-5)) and (u[-1] < 0.999):
        print('GC is implemented only for fully saturated filter.')
        for (name, value) in zip(('L', 's2', 'u_last', 'fl2'),
                                 (L, s2, u[-1], fl2)):
            print('{:6} = {:8g}'.format(name, value))
        print('u = ', u)
        raise ValueError

    y  = model.y
    dy = model.dy
    ds = s2 - s1

    # sample = soil + filter2
    r0 = 0.0

    gc_unsat = (model.porosity * 1/2 * model.density * ds
                * ((r0 + s1)*dy[0]*u[0]
                   + (r0 + s2)*dy[-1]*u[-1]
                   + np.sum((dy[:-1]+dy[1:])
                            *(r0 + s1 + ds*y[1:-1])*u[1:-1])))
    gc_sat   = (1/2 * model.density
                * (model.porosity * (np.power(r0 + s1, 2) - np.power(r0, 2))
                   + (np.power(r0, 2) - np.power(r0 - mass_in, 2))))
    if fl2 > 0.0:
        gc_sat += 1/2 * model.density * model.fp2 * fl2 * (2*(r0 + L) + fl2)

    gc_from_left  = (gc_unsat + gc_sat) / WM_in_tube
    gc_from_right = model.fl1 + L + fl2 - gc_from_left

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
