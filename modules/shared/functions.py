import numpy as np
from math import sqrt

def rpm2radps(x):
    """
      Converts rpm to rad.s^{-1}
    """
    # rpm->rad.s-1:  omega_radps = (2pi)*omega_rps/60
    return x * np.pi/ 30.0

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

def find_omega2g(model, t_total):
    if model.include_acceleration:
        t = t_total - model.t0

        if (model.omega == model.omega_start) or (t > 21.0):
            omega = model.omega
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

            omega_start = model.omega_start
            omega_rel = model.omega - omega_start
            omega = omega_start + omega/omega_base * (model.omega - omega_start)
    else:
        omega = model.omega

    return omega * omega / model.g

def find_omega2g_dec(model, t):
    duration = model.duration
    # omega_end = 0.0, t_end == duration, t in [0, duration]
    omega = (duration - t) / model.deceleration_duration * model.omega

    return omega * omega / model.g

def y2x(y, s1, s2):
    s1_len = np.alen(s1)
    if s1_len != np.alen(s2):
        print('Interfaces array ''s1'' and ''s2'' have to be of the same'
              'lenght. Cannot proceed.')
        exit(1)
    x = np.empty([s1_len, len(y)], float)

    ds = s2 - s1

    for i in range(s1_len):
        x[i, :] = s1[i] + y * ds[i]

    return x

def show_results(extract_data, model, inv_params=None, cov=None):
    from modules.shared.show import ResultsData, DPlots

    data = ResultsData()
    data.extract(extract_data, model, model.params_ref,
                 measurements=model.measurements)
    data.add_value(inv_params=inv_params, cov=cov)

    if model.save_data:
        data.dump(model.experiment_info)

    if model.show_figures:
        dplots = DPlots(data, model.experiment_info)
        dplots.display()

def has_data(x):
    if x is None:
        return False
    elif isinstance(x, np.ndarray):
        return not (x.size == 0)
    else:
        return bool(x)
