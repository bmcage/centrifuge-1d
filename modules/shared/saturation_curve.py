from __future__ import division

"""
  This modules contains function based on Van Genuchten formulas
  and formulas derived from them.
"""
import numpy as np
















P_DEFAULT = np.arange(0, 10000000, 100)

def retention_curve(n, gamma, theta_s, rho, g,
                        theta_r=0.0, p=None, h=None, find_p=True):
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
    if (p is None) and (h is None):
        p = P_DEFAULT

    if h is None:
        h = -10.0* p /rho / g
    elif find_p:
        p = - h * g * rho / 10.0
    else:
        p = None

    theta = theta_r + ((theta_s - theta_r)
                       * h2u(h, n, 1.-1./n, gamma))

    return (p, theta)
