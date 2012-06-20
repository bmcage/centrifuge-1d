import numpy as np

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
    return np.power(10, np.floor(np.log10(np.max(np.abs(v)))))

def scale_array(v, c_coef, result=None):
    """
      Scale the values of array 'v' so that all elements are uniformely divided
      by a scaling factor of 'c_coef'.
      If 'result' is not None, store the resulting values there.
    """
    if c_coef == 1.0: return v

    if result is None:
        return v / c_coef
    else:
        result[:] = v / c_coef
        return result

def f1(t):
    return 1.7032046506 * np.power(t, 1.233644749)
def f2(t):
    return 0.630314472 * np.log(t) + 8.4248850255
def f3(t):
    return 0.1332308098 * np.log(t) + 9.5952480661

def find_omega2g(t_current, model, t_base = 0.0):
    """
    Model includes the acceleration and deceleration of the centrifuge.
    The acceleration model is based on data measured for the centrifuge
    that accelerates to the 600 rpm (i.e. 10 rps) - the evolution of rotational
    speed for other end-speeds is the done by scaling. Deceleration is
    considered to be linear.
    """
    t         = t_current - t_base
    t_end     = model.duration
    omega_max = model.omega
    #print('omg2g: ', model.omega, model._omega)

    if model.include_acceleration:
        if t > t_end:
            if model.deceleration_duration > 0.:
                if t > t_end + model.deceleration_duration:
                    omega = model.omega_end
                else:
                    omega = (model.omega_end
                             + (t_end + model.deceleration_duration - t)
                               / model.deceleration_duration
                               * (omega_max - model.omega_end))
            else:
                omega = omega_max
        elif t > 21.0:
            omega = omega_max
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

            omega = omega/omega_base * omega_max
    else:
        omega = omega_max

        #print('omega: ', omega)
    #print('t = ', t, 'omega = ', omega)
    return np.power(omega, 2)/model.g
    # Previous exponential acceleration:
    #     return (np.power(model.omega_start  + (model.omega - model.omega_start)
    #                      *(1 - np.exp(-model.omega_gamma*t)), 2)
    #             / model.g)

