# -*- coding: utf-8 -*-
"""
Created on Fri Feb 13 14:49:34 2015

@author: benny
"""
from __future__ import print_function, division

import numpy as np

#rpms = [1,500]
#rpms = [500,800,1000,1200]
rpms = [100,500,1000,1300,1650,2000]
#duration = [43200,43200,43200,43200]
#duration = [1800,1800,1800,1800]  #30min
#duration = [120,1800]
duration = [100,78200,14100,8100,14100,8400]
#acceleration_durations = [25,25,25,25]
acceleration_durations = [25,25,25,25,25,25]
#rpms = [500,800,1000,1200,1500,1800,2100,2400,2600,2800,3000]
#duration = [1000,500,500,500,500,500,500,500,500,500,500]
#acceleration_durations = [25,25,25,25,25,30,30,30,30,35,35]

def f1(t):
    """ Helper function for estimating the centrifuge acceleration curve. """
    return 1.7032046506 * np.power(t, 1.233644749)
def f2(t):
    """ Helper function for estimating the centrifuge acceleration curve. """
    return 0.630314472 * np.log(t) + 8.4248850255
def f3(t):
    """ Helper function for estimating the centrifuge acceleration curve. """
    return 0.1332308098 * np.log(t) + 9.5952480661

#we generate rpms output in rad/s
with open("testout_rpms_to2000_SKW.csv", 'w') as f:
    omega_start = 0.;
    timenow = 0.
    for rpm, dur, acceleration_duration in zip(rpms, duration,
                                               acceleration_durations):
        omega_final = rpm *2*np.pi / 60.
        for time in range(dur):
            # Transform t so that acc is in <0, acceleration_duration>
            t = time * 21/acceleration_duration

            if (omega_final == omega_start) or (t > 21.0):
                omega = omega_final
            else:
                omega_base = 10.

                if t > 20.:
                    omega = (21.-t)*f3(t) + (t - 20.)*omega_base
                elif t > 12.:
                    omega = f3(t)
                elif t > 10.:
                    omega = (12.-t)/2.*f2(t) + (t - 10.)/2.*f3(t)
                elif t > 4.:
                    omega =  f2(t)
                elif t > 3.:
                    omega = (4.-t)*f1(t) + (t-3.)*f2(t)
                else:
                    omega = f1(t)

                omega = omega_start + omega/omega_base * (omega_final - omega_start)

            f.write(str(time+timenow) + ';' + str(omega) + '\n')
        timenow += dur
        omega_start = omega_final
