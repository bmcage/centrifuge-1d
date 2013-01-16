from __future__ import division

"""
  This modules contains function based on Van Genuchten formulas
  and formulas derived from them.
"""
import numpy as np

def h2Kh(h, n, m, gamma, Ks, Kh = None):
    tmp1 = np.power(gamma * h, n-1.)
    tmp2 = np.power(1 + gamma * h * tmp1, m/2.)

    if not Kh is None:
        Kh[:] = Ks/tmp2 * np.power(1-tmp1/np.power(tmp2, 2), 2)
    else:
        Kh    = Ks/tmp2 * np.power(1-tmp1/np.power(tmp2, 2), 2)

    return Kh

def u2Ku(u, m, gamma, Ks, Ku = None):
    if not Ku is None:
        Ku[:] = (Ks * np.power(u, 0.5)
                 * np.power(1 - np.power(1 - np.power(u, 1./m), m), 2))
    else:
        Ku    = (Ks * np.power(u, 0.5)
                 * np.power(1 - np.power(1 - np.power(u, 1./m), m), 2))

    return Ku

def dudh(h, n, m, gamma, dudh = None):
    tmp1 = np.power(gamma*h, n-1)
    tmp2 = np.power(1 + gamma*h * tmp1, m+1)

    if not dudh is None:
        dudh[:] = - gamma*(n-1) * tmp1 / tmp2
    else:
        dudh    = - gamma*(n-1) * tmp1 / tmp2

    return dudh

def dhdu(u, n, m, gamma, dhdu = None):
    tmp1 = np.power(u, 1./m)
    tmp2 = np.power(1 - tmp1, m)

    if not dhdu is None:
        dhdu[:] = - 1/gamma/(n-1) / tmp1 / tmp2
    else:
        dhdu    = - 1./gamma/(n-1) / tmp1 / tmp2

    return dudh

def h2u(h, n, m, gamma, u = None):
    if not u is None:
        u[:] = 1/np.power(1+ np.power(gamma * h, n), m)
    else:
        u    = 1/np.power(1+ np.power(gamma * h, n), m)

    return u

def u2h(u, n, m, gamma, h = None):
    if not h is None:
        h[:] =  np.power(np.power(u, -1/m) - 1, 1/n) / gamma
    else:
        h    =  np.power(np.power(u, -1/m) - 1, 1/n) / gamma

    return h

P_DEFAULT = np.arange(0, 10000000, 100)

def retention_curve(n, gamma, theta_s, rho, g,
                        theta_r=0.0, p=P_DEFAULT, h=None):
    """
      Determine the retention curve.

      Parameters:
        n, gamma - van Genuchten soil parameters
        theta_s  - maximal saturation; equals to porosity
        theta_r  - residual saturation
        rho      - fluid density
        g        - gravitation constant
        p        - fluid pressure
        h        - pressure head

      Return values:
        p        - fluid pressure
        theta    - saturation corresponding to pressure p
    """

    if ((p is None) and (h is None)) or ((not p is None) and (not h is None)):
        print('To display retention curve either pressure p or hydraulic head '
              'h must be set, but not both.')
        print('h = ', h)
        print('p = ', p)
        return None

    if h is None:
        h = -10.0* p /rho / g
    else:
        p = - h * g * rho / 10.0

    theta = theta_r + ((theta_s - theta_r)
                       * h2u(h, n, 1.-1./n, gamma))

    return (p, theta)
