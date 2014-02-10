from __future__ import division

"""
  This modules contains saturation curve object(s) for single flow
"""
import numpy as np
from functions import MonoCubicInterp

__NR = 400
P_DEFAULT = np.linspace(-1, 9, __NR)
P_DEFAULT = np.power(10* np.ones(__NR), P_DEFAULT)
P_DEFAULT[0] = 0
#P_DEFAULT = np.arange(0, 10000000, 100)

########################################################################
#                             Common class                             #
########################################################################

class SC_base():
    def retention_curve(self, theta_s, rho, g, theta_r=0.0, p=None, h=None,
                        find_p=True):
        """
          Determine the retention curve.

          Parameters:
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

        theta = theta_r + (theta_s - theta_r) * self.h2u(h)

        return (p, theta)

    def retention_curve_u(self, h=None, p=None, g=None, rho=None):
        """ Retention curve in terms of relative saturation u. """

        if (p is None) and (h is None):
            assert (g is not None) and (rho is not None)
            p = P_DEFAULT

        if h is None:
            h = -10.0* p /rho / g

        return (h, self.h2u(h))

    def conductivity_curve(self, Ks, theta_s, theta_r=0.0, u=None, p=None,
                           h=None, rho=None, g=None):
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
            u = self.h2u(h)
        elif not p is None:
            if (rho is None) or (g is None):
                print("Conductivity curve: neither 'rho' nor 'g' can be 'None' !"
                      "No conductivity curve is computed...")
                return ([], [])
            u = self.h2u(-10.0* p /rho / g)

        theta = theta_r + (theta_s - theta_r) * u
        K     = self.u2Ku(u, Ks)

        return (theta, K)

    def conductivity_curve_u(self, Ks, theta_s, theta_r=0.0, u=None, p=None,
                           h=None, rho=None, g=None):
        th, K = self.conductivity_curve(Ks, theta_s, theta_r, u, p,
                                   h, rho, g)
        u = (th-theta_r)/(theta_s - theta_r)
        return (u, K)

    def get_dyn_h_init(self, c_gammah, h_init_max):
        """
        Initial values of h~0 can cause troubles to the solver to start
        depending on parameters of the SC class. To ensure "smooth" start
        we compute a dynamically obtained 'h_init' value based on
        actual values. This may be important in the parameters
        optimization process.
        """
        pass

    def add_transformations_fns(self, transform, untransform, max_value):
        """
        When perfoming optimization problem over the SC parameters,
        we can work directly with the parameters values or
        some transformation of them (e.g. to ensure we keep the
         searching interval bounded).
        This parameters transformation is performed ONLY to determine
        next value in the inverse optimization process, i.e.methods
        using these parameters always obtain untransformed values.
        """
        pass


########################################################################
#                       van Genuchten model                            #
########################################################################

class SC_vanGenuchten(SC_base):
    """ van Genuchten model.
    This is given by writing effective saturation u (=S_e) in terms of
    pressure head h, and relative permeability in terms of S_e using the same
    parameters.

    We have

         S_e = u =  1 / (1 + (\gamma h)^n)^m,   m = 1 - 1/n

         k(u)  =  u^0.5 (1 - (1 - u^{1/m})^m )^2
    """
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
        """
        For a given head h and saturated hydraulic conductivity Ks, what is
        K(h) = Ks k(u(h)).
        If Kh given, it is used to store the result
        """
        tmp1 = np.power(self._gamma * h, self._n-1.)
        tmp2 = np.power(1 + self._gamma * h * tmp1, self._m/2.)

        if not Kh is None:
            Kh[:] = Ks/tmp2 * np.power(1-tmp1/np.power(tmp2, 2), 2)
        else:
            Kh    = Ks/tmp2 * np.power(1-tmp1/np.power(tmp2, 2), 2)

        return Kh

    def u2Ku(self, u, Ks, Ku = None):
        """
        For a given effective saturation u and saturated hydraulic conductivity
        Ks, what is K(u) = Ks k(u)
        """
        m = self._m

        if not Ku is None:
            Ku[:] = (Ks * np.power(u, 0.5)
                     * np.power(1 - np.power(1 - np.power(u, 1./m), m), 2))
        else:
            Ku    = (Ks * np.power(u, 0.5)
                     * np.power(1 - np.power(1 - np.power(u, 1./m), m), 2))

        return Ku

    def dudh(self, h, dudh = None):
        """
        Return value of du/dh in h
        """
        n     = self._n
        gamma = self._gamma
        m = self._m

        tmp1 = np.power(gamma*h, n-1)
        tmp2 = np.power(1 + gamma*h * tmp1, m+1)

        if not dudh is None:
            dudh[:] = - gamma*(n-1) * tmp1 / tmp2
        else:
            dudh    = - gamma*(n-1) * tmp1 / tmp2

        return dudh

    def dhdu(self, u, dhdu = None):
        """
        Return value dh/du in terms of u for the given n, m and gamma
        """
        n     = self._n
        gamma = self._gamma
        m = self._m

        tmp1 = np.power(u, 1./m)
        tmp2 = np.power(1 - tmp1, m)

        if not dhdu is None:
            dhdu[:] = - 1/gamma/(n-1) / tmp1 / tmp2
        else:
            dhdu    = - 1./gamma/(n-1) / tmp1 / tmp2

        return dhdu

    def h2u(self, h, u = None):
        """
        Return the value of effective saturation u corresponding with h.
        """
        if not u is None:
            u[:] = 1/np.power(1+ np.power(self._gamma * h, self._n), self._m)
        else:
            u    = 1/np.power(1+ np.power(self._gamma * h, self._n), self._m)

        return u

    def u2h(self, u, h = None):
        """
        Return the value of h for the given effective saturation u and
        parameters passed
        """
        n     = self._n
        gamma = self._gamma
        m = self._m

        if not h is None:
            h[:] =  (np.power(np.power(u, -1./m) - 1., 1./n)
                     / gamma)
        else:
            h    =  (np.power(np.power(u, -1./m) - 1, 1/n)
                     / gamma)

        return h

    def get_dyn_h_init(self, c_gammah, h_init_max):
        """
        Dynamically set 'h_init' value based on the 'gamma' parameter.
        """
        return  min(c_gammah / self._gamma, h_init_max)

    def add_transformations_fns(self, transform, untransform, max_value):
        """
        Transform/untransform methods for 'n' and 'gamma' parameters
        """
        transform['n']   = lambda n: max(np.log(n - 1.0), -max_value)
        untransform['n'] = lambda n_transf: 1+min(np.exp(n_transf), max_value)

        transform['gamma']   = lambda gamma: max(np.log(-gamma), -max_value)
        untransform['gamma'] = lambda gamma_transf: -min(np.exp(gamma_transf), max_value)

########################################################################
#                           Free-form  model                           #
########################################################################

class SC_freeform(SC_base):
    """
    A freefrom description of the saturation curve.
    In this, a discrete version of h is used: {h_i}, with -inf < h_i < 0.
    For these gridpoints, the value of effective saturation is passed:
        u_i = u(h_i)
    and also the value of the relative permeability:
        k_i = k(h_i)

    We then determine u(h) and k(h) as a monotone cubic interpolation spline
    of the given points {(h_i, u_i)} and {(h_i, k_i)}
    From this cubic spline the derivatives can be extracted
    """
    def __init__(self, hi=None, ui=None, ki=None):
        if hi is not None or ui is not None or ki is not None:
            SC_freeform.check_params(hi, ui, ki)
            self.__set_values(hi, ui, ki)
        else:
            self._hi = None
            self._ui = None
            self._ki = None

    def __set_values(self, hi, ui, ki):
        #are internal parameters are such that all is monotone increasing
        # in terms of h
        self._hi = hi
        self._ui = ui
        self._ki = ki
        self._hnodes = np.empty(len(hi)+2, float)
        self._uvals = np.empty(len(hi)+2, float)
        self._kvals = np.empty(len(hi)+2, float)
        #we use the log values  for head.
        self._hnodes[1:-1] = -1*np.log(-hi[::-1])
        self._hnodes[0] = -1*np.log(-hi[0]*1e3)
        self._hnodes[-1] - -1*np.log(1e-28)
        #with log for head, saturation is linear
        self._uvals[1:-1] = ui[:]
        self._uvals[0] = 0.
        self._uvals[-1] = 1
        #with log for head, also rel perm should be log values
        self._kvals[1:-1] = np.log(ki[:])
        self._kvals[0] = 0.
        self._kvals[-1] = 1

        # we now construct two cubic monotone interpolations: U(H) and K(H)
        self.logh2u = MonoCubicInterp(self._hnodes, self._uvals)
        self.logh2logk = MonoCubicInterp(self._hnodes, self._kvals)

    @staticmethod
    def check_params(hi, ui, ki):
        if not hi or not ui or not ki:
            raise Exception, 'Some parameters are not given: '\
                    'hi:%s, ui:%s, ki:%s' % (str(hi), str(ui), str(ki))
        if not (len(hi)==len(ui)==len(ki)):
            raise Exception, 'Parameters must have equal length'
        #hi must be monotone ascending > 0, and ui, ki monotone decreasing
        ho = hi[0]
        uo = ui[0]
        ko = ki[0]
        for h,u,k in zip(hi[1:],ui[1:],ki[1:]):
            if not h>ho or not h<0.:
                raise Exception, 'Hydraulic head h must be negative and a '\
                    'monotone ascending array, instead %s' % str(hi)
            if not u>uo or not uo>0:
                raise Exception, 'Effective saturation Se must be positive and a '\
                    'monotone ascending array in terms of h, instead %s' % str(ui)
            if not k<ko or not ko>0:
                raise Exception, 'Relative permeability k must be positive and a '\
                    'monotone ascending array in terms of h, instead %s' % str(ki)

        #The edges of ui and ki are fixed and will not be optimized
        if not (ui[0] >0) or not (ui[-1] < 1):
            raise Exception, 'Effective saturation Se starts at 0, ends at 1, '\
                    'only pass values between these extremes!'
        if not (ki[0] > 0) or not (ki[-1] < 1):
            raise Exception, 'Relative permeability k starts at 0, ends at 1, '\
                    'only pass values between these extremes!'

    def set_parameters(self, params):
        """
        Setting the parameters. Only ui and ki are parameters!
        """
        if 'ui' in params:
            ui = params['ui']
        else:
            ui = self._ui
        if 'ki' in params:
            ki = params['ki']
        else:
            ki = self._ki
        SC_freeform.check_params(self._hi, ui, ki)
        self.__set_values(self._hi, ui, ki)

    def h2Kh(self, h, Ks, Kh = None):
        """
        For a given head h and saturated hydraulic conductivity Ks, what is
        K(h) = Ks k(u(h)).
        If Kh given, it is used to store the result
        """
        tmp1 = np.exp( self.logh2logk(-np.log(-h) ))

        if not Kh is None:
            Kh[:] = Ks * tmp1
        else:
            Kh    = Ks * tmp1
        return Kh

    def u2Ku(self, u, Ks, Ku = None):
        """
        For a given effective saturation u and saturated hydraulic conductivity
        Ks, what is K(u) = Ks k(u)
        """
        h = self.u2h(u)
        return self.h2Kh(h, Ks, Ku)

    def h2u(self, h, u = None):
        """
        Return the value of effective saturation u corresponding with h.
        """
        tmp1 = np.exp( self.logh2u(-np.log(-h) ))
        if not u is None:
            u[:] = tmp1[:]
        else:
            u    = tmp1

        return u

    def u2h(self, u, h = None):
        """
        Return the value of h for the given effective saturation u and
        parameters passed
        """
        tmph = self.logh2u.root(u)
        tmph = - np.exp(-tmph)

        if not h is None:
            h[:] =  tmph[:]
        else:
            h    =  tmph

        return h

    def dudh(self, h, dudh = None):
        """
        Return value of du/dh in h
        We have logh2u which is u(logh(h)) so
        du/dh = d logh2u / dlogh .  dlogh/dh, with dlogh/dh = -1/h
        """
        logh = -np.log(-h)
        tmp1 = self.logh2u.derivative(logh,der=1) * -1/h

        if not dudh is None:
            dudh[:] = tmp1[:]
        else:
            dudh    = tmp1

        return dudh
