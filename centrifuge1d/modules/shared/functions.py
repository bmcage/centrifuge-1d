from __future__ import division, print_function

import numpy as np
from math import sqrt
from sys import stdout

def rpm2radps(omega):
    """
      Converts rpm to rad/s
    """
    # rpm->rad.s-1:  omega_radps = (2pi)*omega_rpm/60
    _rpm2radps = lambda omega: omega * np.pi/ 30.0

    if np.isscalar(omega) or type(omega) is np.ndarray:
        omega_converted = _rpm2radps(omega)
    elif type(omega) in (list, tuple):
        omega_converted = [_rpm2radps(omg) for omg in omega]
    else:
        raise ValueError('Unknown type for omega')

    return omega_converted

def radps2rpm(omega):
    """
      Converts rad/s to rpm
    """
    # rad/s->rpm:  omega_rpm = 60*omega_radps/(2pi)
    _radps2rpm = lambda omega: omega * 30/np.pi

    if np.isscalar(omega) or type(omega) is np.ndarray:
        omega_converted = _radps2rpm(omega)
    elif type(omega) in (list, tuple):
        omega_converted = [_radps2rpm(omg) for omg in omega]
    else:
        raise ValueError('Unknown type for omega')

    return omega_converted

def lagrangian_derivative_coefs(dx):
    """
    Returns the coeficients for the Lagrangeand derivative of the differences
    array 'dx'. The first point has a right derivative, last point a left
    derivative and central difference is for the mid-points.
    Typical usage: res the values over grid x, then
     dresdx[0]    = ldc1[0]   *  res[0] + ldc2[0]   *   res[1] + ldc3[0]  * res[2]
     dresdx[1:-1] = ldc1[1:-1]*res[:-2] + ldc2[1:-1]*res[1:-1] + ldc3[1:-1]*res[2:]
     dresdx[-1]   = ldc1[-1]  * res[-3] + ldc2[-1]  *  res[-2] + ldc3[-1] * res[-1]
    """
    ldc1 = np.concatenate(([-(2*dx[0]+dx[1])/(dx[0]*(dx[0]+dx[1]))],
                          -dx[1:]/(dx[:-1]*(dx[:-1]+dx[1:])),
                          [dx[-1]/(dx[-2]*(dx[-2]+dx[-1]))]))
    ldc2 = np.concatenate(([(dx[0]+dx[1])/(dx[1]*dx[0])],
                          (dx[1:] - dx[:-1])/dx[:-1]/dx[1:],
                          [-(dx[-1]+dx[-2])/(dx[-2]*dx[-1])]))
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

def f1(t):
    """ Helper function for estimating the centrifuge acceleration curve. """
    return 1.7032046506 * np.power(t, 1.233644749)
def f2(t):
    """ Helper function for estimating the centrifuge acceleration curve. """
    return 0.630314472 * np.log(t) + 8.4248850255
def f3(t):
    """ Helper function for estimating the centrifuge acceleration curve. """
    return 0.1332308098 * np.log(t) + 9.5952480661

def find_omega2g_fh(r0_fall):
    """
      Function for determination of the omega^2/g coef when simulation
      falling-head test (i.e. simulating 1g environment).
    """
    return 1./r0_fall

def find_omega2g(t, omega_final, omega_start, g,
                 include_acceleration, acceleration_duration):
    """
      Function for determinantion of the omega^2/g coef under centrifugation.
      This includes also the acceleration curve (but not deceleration curve).
      Time 't' is the time since the centrifugation (phase) started.
    """
    if include_acceleration:
        # Transform t so that acc is in <0, acceleration_duration>
        t = t * 21/acceleration_duration

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

            omega_rel = omega_final - omega_start
            omega = omega_start + omega/omega_base * (omega_final - omega_start)
    else:
        omega = omega_final

    return omega * omega / g

def find_omega2g_dec(t, deceleration_duration, omega_start, g):
    """
      Function for determinantion of the omega^2/g coef (under centrifugation)
      when decelerating.
    """
    # omega_end = 0.0, t_end == duration, t in [0, duration]
    omega = ((deceleration_duration - t)
             / deceleration_duration * omega_start)

    return omega * omega / g

def y2x(y, s1, s2):
    """ Transform interval y=<0,1> to original interval x. """
    s1_len = np.alen(s1)
    if s1_len != np.alen(s2):
        print('Interfaces array ''s1'' and ''s2'' have to be of the same'
              'lenght. Cannot proceed.')
        exit(1)
    x = np.empty([s1_len, len(y)], float)

    ds = s2 - s1

    if s1_len > 1:
        for i in range(s1_len):
            x[i, :] = s1[i] + y * ds[i]
    else:
        x[:] = s1 + y*ds

    return x

def has_data(x):
    if x is None:
        return False
    elif isinstance(x, np.ndarray):
        return not (x.size == 0)
    else:
        return bool(x)

def phases_end_times(a_duration, d_duration, g_duration,
                     include_acceleration):
    """
      Calculate the stop-times for the direct solver based on the information
      of the individual phases duration
      (i.e. acceleration/centrifugation + deceleration + sample left under 1g).
    """

    if a_duration is None:
        a_duration = 0.0
    else:
        a_duration = np.asarray(a_duration, dtype=float)

    if include_acceleration and (not d_duration  is None):
        a_duration += d_duration

    if g_duration is None:
        g_duration = 0.0
    else:
        g_duration = np.asarray(g_duration, dtype=float)

    duration_times = a_duration + g_duration
    if not np.any(duration_times): return None # no times were specified

    if np.isscalar(duration_times):
        duration_times = duration_times.reshape([1,])

    stop_times = np.cumsum(duration_times)

    return stop_times

def compare_data(name, value_computed, value_measured = None,
                 stream=None):
    """
      Display the measured and computed values, abs, rel and LSQ error.
      Data is written to 'stream', which is by default stdout.
    """

    if stream is None: stream=stdout

    name_len = len(name)
    data_computed = np.asarray(value_computed, dtype=float)
    disp_all = (not value_measured is None)

    if disp_all:
        data_measured = np.asarray(value_measured, dtype=float)

        measured_filter = (data_measured == 0.0)
        norm_measured = np.abs(data_measured)
        norm_measured[measured_filter] = 1e-60
        rel_error = (data_computed - data_measured) / norm_measured * 100.
        rel_error[abs(rel_error) > 1e50] = np.inf
        abs_error = np.abs(data_computed - data_measured)

    i0 = 0
    in_row = 10
    remaining = np.alen(data_computed)

    float_disp_size = 12
    fstr = '% {}.6f'.format(float_disp_size)

    print('\n', file=stream)
    while remaining > 0:
        if remaining > in_row:
            disp_items = in_row
        else:
            disp_items = remaining

        print('%s measured: ' % name,
              disp_items * fstr % tuple(data_measured[i0:i0+disp_items]),
              file=stream)
        if disp_all:
            print('%s computed: ' % name,
                  disp_items * fstr % tuple(data_computed[i0:i0+disp_items]),
                  file=stream)
            print('AbsError: ', name_len * ' ',
                  disp_items * fstr % tuple(abs_error[i0:i0+disp_items]),
                  file=stream)
            print('Error (%):', name_len * ' ',
                  disp_items * fstr % tuple(rel_error[i0:i0+disp_items]),
                  file=stream)

        remaining = remaining - disp_items
        print((16 + float_disp_size*in_row) * '-', file=stream)
        i0 = i0 + in_row

    LSQ_error = np.sum(np.power(data_computed - data_measured, 2))
    print('LSQ error ' + "'" + name + "':", LSQ_error, file=stream)

    return LSQ_error

#*******************************************************
#
#  Averaging functions
#     From: http://www.swharden.com/blog/2008-11-17-linear-data-smoothing-in-python/
#     From: http://pandas.pydata.org
#*******************************************************

def smoothing_linear(vlist, degree=5):
    # smooth based on linear averaging before after of window=2*degree-1 points
    window=degree*2-1
    smoothed = np.zeros(len(vlist), float)
    smoothed[:degree-1] = vlist[:degree-1]
    for i in range(len(vlist)-window):
        smoothed[degree-1+i]=sum(vlist[i:i+window])/window
    smoothed[-degree:] = vlist[-degree:]
    return smoothed

def smoothing_triangle(vlist, degree=5):
    # smooth based on triangle averaging before after of window=2*degree-1 points
    weight=[]
    window=degree*2-1
    smoothed = np.zeros(len(vlist), float)
    smoothed[:degree-1] = vlist[:degree-1]
    for x in range(1, 2*degree):weight.append(degree-abs(degree-x))
    w=np.array(weight)
    smoothed[:degree-1] = vlist[:degree-1]
    for i in range(len(vlist)-window):
        smoothed[degree-1+i]=sum(np.array(vlist[i:i+window])*w)/sum(w)
    smoothed[-degree:] = vlist[-degree:]
    return smoothed

def smoothing_gaussian(vlist, degree=5):
    # smooth based on gaussian averaging before after of window=2*degree-1 points
    window=degree*2-1
    weight=np.array([1.0]*window)
    weightGauss=[]
    for i in range(window):
        i=i-degree+1
        frac=i/float(window)
        gauss=1/(np.exp((4*(frac))**2))
        weightGauss.append(gauss)
    weight=np.array(weightGauss)*weight
    smoothed = np.zeros(len(vlist), float)
    smoothed[:degree-1] = vlist[:degree-1]
    for i in range(len(vlist)-window):
        smoothed[degree-1+i]=sum(np.array(vlist[i:i+window])*weight)/sum(weight)
    smoothed[-degree:] = vlist[-degree:]
    return smoothed

def smoothing_gaussian_rec(vlist, rec = 2, degree=5):
    # smooth based on recursive gaussian averaging, averaging starts before and
    # after of window=2*degree-1 points
    assert rec > 1
    for entry in range(rec):
        vlist = smoothing_gaussian(vlist, degree)
    return vlist

def smoothing_gaussian_exp(vlist, degree=5, smooth_fac=0.5):
    #exponential smoothing based on underlying gaussian smoothing
    #smooth_fac is the exponential smoothing factor. If 0 no extra smoothing, 
    assert 0 < smooth_fac < 1
    vlist = smoothing_gaussian(vlist, degree)
    smoothed = np.empty(len(vlist), float)
    smoothed[:] = vlist[:]
    #now exponential smoothing on this
    for i in range(1, len(vlist)):
        smoothed[i]=smooth_fac * smoothed[i-1] + (1-smooth_fac) * smoothed[i]
    return smoothed

def smoothing_exp(vlist, smooth_fac=0.5):
    #exponential smoothin
    #smooth_fac is the exponential smoothing factor. If 0 no extra smoothing
    assert 0 < smooth_fac < 1
    smoothed = np.empty(len(vlist), float)
    smoothed[:] = vlist[:]
    #now exponential smoothing on this
    for i in range(1, len(vlist)):
        smoothed[i]=smooth_fac * smoothed[i-1] + (1-smooth_fac) * smoothed[i]
    return smoothed
