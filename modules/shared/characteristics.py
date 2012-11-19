from __future__ import division

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

def calc_gc(u, y, dy, s1, s2, mass_in, rs_sat, re_sat, soil_porosity, fl2, fp2,
            fr2, WM_in_tube, fluid_density, from_end=None):
    """
      Determine the gravitational center of water in the sample. The water
      on the inflow is taken into account (but not water on the outflow).
      Computed GC is measured from the BEGINNING of the soil sample
      (in the direction from centrifuge axis)

      Arguments:
      u - relative saturation
      y - transformed interval where holds: u(x) = u(fl1 + s1 + (s2-s1)y)
      dy - difference of y: dy = diff(y) = y(i+1) - y(i)
      mass_in - water in the inflow chamber
      rs_sat - starting radius of saturated part
      re_sat - ending radius of saturated part
      s1, s2 - interfaces (s1 < s2)
      fl2, fp2, fr2 - ending filter length, porosity and distance from sample
                      beginning (to filter's beginning)
      WM_in_tube - amount of water contained inside the tube
      from_end - if specified, computed GC will be returned as distance
                 from the "from_end" point
    """

    ds  = s2 - s1
    r0 = 0.0

    gc_unsat = (soil_porosity * 1/2 * fluid_density * ds
                * ((r0 + s1)*dy[0]*u[0]
                   + (r0 + s2)*dy[-1]*u[-1]
                   + np.sum((dy[:-1]+dy[1:])
                            *(r0 + s1 + ds*y[1:-1])*u[1:-1])))
    gc_sat   = \
      (1/2 * fluid_density
       * (soil_porosity * (np.power(r0 + re_sat, 2) - np.power(r0 + rs_sat, 2))
          + (np.power(r0, 2) - np.power(r0 - mass_in, 2))))
    if fl2 > 0.0:
        gc_sat += 1/2 * fluid_density * fp2 * fl2 * (2*(r0 + fr2) + fl2)

    gc = (gc_unsat + gc_sat) / WM_in_tube

    if not from_end is None: gc = from_end - gc

    return gc

def calc_rm(t, u, mass_in, mass_out, s1, s2, model):

    raise NotImplementedError('Calculation of rotational momentum is not '
                              'verified and therefore not provided.')

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
