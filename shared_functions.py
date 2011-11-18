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
    ldc2 = np.concatenate(([dx[0]+dx[1]/(dx[1]*dx[0])],
                          (dx[1:] - dx[:-1])/dx[:-1]/dx[1:],
                          [(dx[-1]+dx[-2])/(dx[-2]*dx[-1])]))
    ldc3 = np.concatenate(([-dx[0]/(dx[1]*(dx[1]+dx[0]))],
                           dx[:-1]/(dx[1:]*(dx[:-1]+dx[1:])),
                           [(2*dx[-1]+dx[-2])/(dx[-1]*(dx[-2]+dx[-1]))]))

    return ldc1, ldc2, ldc3

def h2Kh(h, n, m, gamma, Ks, Kh = None, tmp1 = None, tmp2 = None):
    if dudh and tmp1 and tmp2:
        tmp1[:] = np.power( gamma * h, n-1.)
        tmp2[:] = np.power(1 + gamma * h * tmp1, m/2.)
        Kh[:]   = Ks/tmp2 * np.power(1-tmp1/np.power(tmp2, 2), 2)
    else:
        tmp1  = np.power( gamma * h, n-1.)
        tmp2 = np.power(1 + gamma * h * tmp1, m/2.)
        Kh   = Ks/tmp2 * np.power(1-tmp1/np.power(tmp2, 2), 2)

    return Kh

def dudh(h, n, m, gamma, Ks, dudh = None, tmp1 = None, tmp2 = None):
    if dudh and tmp1 and tmp2:
        tmp1[:] = np.power(gamma*h, n-1)
        tmp2[:] = np.power(1 + gamma*h * tmp1, m+1)
        dudh[:] = - gamma*(n-1) * tmp1 / tmp2
    else:
        tmp1 = np.power(gamma*h, n-1)
        tmp2 = np.power(1 + gamma*h * tmp1, m+1)
        dudh = - gamma*(n-1) * tmp1 / tmp2

    return dudh

def h2u(h, n, m, gamma, u = None):
    if u:
        u[:] = 1/np.power(1+ np.power(gamma * h, n), m)
    else:
        u    = 1/np.power(1+ np.power(gamma * h, n), m)

    return u

def u2h(u, n, m, gamma, h = None):
    if h:
        h[:] =  np.power(np.power(u, -1/m) - 1, 1/n) / gamma
    else:
        h    =  np.power(np.power(u, -1/m) - 1, 1/n) / gamma

    return h
