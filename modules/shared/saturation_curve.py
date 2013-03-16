from __future__ import division

"""
  This modules contains saturation curve object(s)
"""
import numpy as np

########################################################################
#                       van Genuchten model                            #
########################################################################

class SC_vanGenuchten():
    def __init__(self, n=None, gamma=None):
        self._n = n
        if n is None:
            self._m = None
        else:
            self._m = 1. - 1./n
        self._gamma = gamma

    def set_parameters(self, params):
        if 'n' in params:
            n = params['n']
            self._n = n
            if n is None:
                self._m = None
            else:
                self._m = 1. - 1./n

        if 'gamma' in params:
            self._gamma = params['gamma']

    def h2Kh(self, h, Ks, Kh = None):
        tmp1 = np.power(self._gamma * h, self._n-1.)
        tmp2 = np.power(1 + self._gamma * h * tmp1, self._m/2.)

        if not Kh is None:
            Kh[:] = Ks/tmp2 * np.power(1-tmp1/np.power(tmp2, 2), 2)
        else:
            Kh    = Ks/tmp2 * np.power(1-tmp1/np.power(tmp2, 2), 2)

        return Kh

    def u2Ku(self, u, Ks, Ku = None):
        m = self._m

        if not Ku is None:
            Ku[:] = (Ks * np.power(u, 0.5)
                     * np.power(1 - np.power(1 - np.power(u, 1./m), m), 2))
        else:
            Ku    = (Ks * np.power(u, 0.5)
                     * np.power(1 - np.power(1 - np.power(u, 1./m), m), 2))

        return Ku

    def dudh(self, h, dudh = None):
        n     = self._n
        gamma = self._gamma

        tmp1 = np.power(gamma*h, n-1)
        tmp2 = np.power(1 + gamma*h * tmp1, self._m+1)

        if not dudh is None:
            dudh[:] = - gamma*(n-1) * tmp1 / tmp2
        else:
            dudh    = - gamma*(n-1) * tmp1 / tmp2

        return dudh

    def dhdu(self, u, n, m, gamma, dhdu = None):
        tmp1 = np.power(u, 1./m)
        tmp2 = np.power(1 - tmp1, m)

        if not dhdu is None:
            dhdu[:] = - 1/gamma/(n-1) / tmp1 / tmp2
        else:
            dhdu    = - 1./gamma/(n-1) / tmp1 / tmp2

        return dudh

    def h2u(self, h, u = None):
        if not u is None:
            u[:] = 1/np.power(1+ np.power(self._gamma * h, self._n), self._m)
        else:
            u    = 1/np.power(1+ np.power(self._gamma * h, self._n), self._m)

        return u

    def u2h(self, u, n, m, gamma, h = None):
        if not h is None:
            h[:] =  (np.power(np.power(u, -1./self._m) - 1., 1./self._n)
                     / self._gamma)
        else:
            h    =  (np.power(np.power(u, -1./self._m) - 1, 1/self._n)
                     / self._gamma)

        return h

    def get_dyn_h_init(self, c_gammah, h_init_max):
        return  min(c_gammah / self._gamma, h_init_max)

    def add_transformations_fns(self, transform, untransform, max_value):
        transform['n']   = lambda n: max(np.log(n - 1.0), -max_value)
        untransform['n'] = lambda n_transf: 1+min(np.exp(n_transf), max_value)

        transform['gamma']   = lambda gamma: max(np.log(-gamma), -max_value)
        untransform['gamma'] = lambda gamma_transf: -min(np.exp(gamma_transf), max_value)

########################################################################
#                           Free-form  model                           #
########################################################################



########################################################################
#                           Common utilities                           #
########################################################################

__NR = 400
P_DEFAULT = np.linspace(-1, 9, __NR)
P_DEFAULT = np.power(10* np.ones(__NR), P_DEFAULT)
P_DEFAULT[0] = 0
#P_DEFAULT = np.arange(0, 10000000, 100)

def retention_curve(SC, theta_s, rho, g, theta_r=0.0, p=None, h=None,
                    find_p=True):
    """
      Determine the retention curve.

      Parameters:
        SC       - saturation curve object
        theta_s  - maximal saturation; equals to porosity
        theta_r  - residual saturation
        rho      - fluid density
        g        - gravitation constant
        p        - fluid pressure
        h        - pressure head
        find_p   - if True and h is supplied, the p is determined
                   (otherwise p=None is returned)

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

    theta = theta_r + (theta_s - theta_r) * SC.h2u(h)

    return (p, theta)

def conductivity_curve(SC, Ks, theta_s, theta_r=0.0, u=None, p=None, h=None,
                       rho=None, g=None):
    """
      Determine the conductivity curve.

      Parameters:
        SC       - saturation curve object
        Ks       - saturated hydraulic conductivity
        theta_s  - maximal saturation; equals to porosity
        theta_r  - residual saturation
        rho      - fluid density
        g        - gravitation constant
        p        - fluid pressure
        h        - pressure head
        u        - relative saturation

      Return values:
        theta    - saturation corresponding to pressure p/relative saturation u/
                   hydraulic head h
        K        - conductivity
    """

    if (u is None) and (p is None) and (h is None):
        p = P_DEFAULT

    if not h is None:
        u = SC.h2u(h)
    elif not p is None:
        if (rho is None) or (g is None):
            print("Conductivity curve: neither 'rho' nor 'g' can be 'None' !"
                  "No conductivity curve is computed...")
            return ([], [])
        u = SC.h2u(-10.0* p /rho / g)

    theta = theta_r + (theta_s - theta_r) * u
    K     = SC.u2Ku(u, Ks)
    print('KKKKK', Ks, K)

    return (theta, K)
