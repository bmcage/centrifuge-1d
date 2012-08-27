import numpy as np
from math import sqrt

def rpm2radps(x):
    """
      Converts rpm to rad.s^{-1}
    """
    # rpm->rad.s-1:  omega_radps = (2pi)*omega_rps/60
    return x * np.pi/ 30.0

def measurements_time(model):
    t_duration = model.get_iterable_value('duration')
    if not t_duration:
        t_duration = model.duration
    t_duration = np.asarray(t_duration, dtype=float)

    if model.include_acceleration:
        t_duration[:] = t_duration + model.deceleration_duration

    t_fh_duration = model.get_iterable_value('fh_duration')
    if not t_fh_duration:
        t_fh_duration = model.fh_duration

    return np.cumsum(t_duration + np.asarray(t_fh_duration, dtype=float))

def lagrangean_derivative_coefs(dx):
    """
    Returns the coeficients for the Lagrangeand derivative of the differences
    array 'dx'. The first point has a right derivative, last point a left
    derivative and central difference is for the mid-points.
    """
    ldc1 = np.concatenate(([-(2*dx[0]+dx[1])/(dx[0]*(dx[0]+dx[1]))],
                          -dx[1:]/(dx[:-1]*(dx[:-1]+dx[1:])),
                          [dx[-1]/(dx[-2]*(dx[-2]+dx[-1]))]))
    ldc2 = -np.concatenate(([(dx[0]+dx[1])/(dx[1]*dx[0])],
                          (dx[1:] - dx[:-1])/dx[:-1]/dx[1:],
                          [(dx[-1]+dx[-2])/(dx[-2]*dx[-1])]))
    ldc3 = np.concatenate(([-dx[0]/(dx[1]*(dx[1]+dx[0]))],
                           dx[:-1]/(dx[1:]*(dx[:-1]+dx[1:])),
                           [(2*dx[-1]+dx[-2])/(dx[-1]*(dx[-2]+dx[-1]))]))

    return ldc1, ldc2, ldc3

def right_derivative(dx12, fx13):
    [dx1, dx2]      = dx12
    [fx1, fx2, fx3] = fx13

    derivative = (dx2/(dx1*(dx1+dx2)) * fx1
                  - (dx2 + dx1)/(dx1 * dx2) * fx2
                  + (2*dx2 + dx1)/(dx2*(dx1+dx2)) * fx3)

    return derivative

def determine_scaling_factor(v):
    return np.power(10, -np.floor(np.log10(np.max(np.abs(v)))))

def f1(t):
    return 1.7032046506 * np.power(t, 1.233644749)
def f2(t):
    return 0.630314472 * np.log(t) + 8.4248850255
def f3(t):
    return 0.1332308098 * np.log(t) + 9.5952480661

def find_omega2g_fh(model, t):
    return 1/model.r0_fall

def find_omega2g(model, t):
    if model.include_acceleration:
        if t > 21.0:
            omega = model.omega_max
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

            omega = omega/omega_base * model.omega_max
    else:
        omega = model.omega_max

    return omega * omega / model.g

def find_omega2g_dec(model, t):
    duration = model.duration
    # omega_end = 0.0, t_end == duration, t in [0, duration]
    omega = (duration - t) / deceleration_duration * model.omega_max)

    return omega * omega / model.g
