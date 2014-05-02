from __future__ import division, print_function

"""
  This modules contains saturation curve object(s) for single flow
"""
import numpy as np
from .interpolate import MonoCubicInterp, QuadraticBspline, PiecewiseLinear

__NR = 400
P_DEFAULT = np.linspace(-1, 9, __NR)
P_DEFAULT = np.power(10* np.ones(__NR), P_DEFAULT)
P_DEFAULT[0] = 0
#P_DEFAULT = np.arange(0, 10000000, 100)

SC_vG     = 1
SC_FF_CUB = 2
SC_FF_BS  = 3
SC_FF_LIN = 4

def get_dict_value(dictionary, keys):
    if type(keys) is str:                 # single item
        return dictionary[keys]
    else:                                 # list of keys
        return [dictionary[key] for key in keys]

def create_SC(data):
    """
    Returns a new SC object. Data is expected to be either a dict with values
    or an instance of Configuration class.
    """
    if type(data) is dict:                # dict
        get_value = lambda keys: get_dict_value(data, keys)
    else:                                 # Configuration class instance
        get_value = data.get_value

    SC_type = get_value('sc_type')
    if SC_type == SC_vG:
        SC = SC_vanGenuchten(get_value('n'), get_value('gamma'))
    elif SC_type in (SC_FF_CUB, SC_FF_BS, SC_FF_LIN):
        (hi, ui, ki, max_refine) = get_value(('hi', 'ui', 'ki', 'sc_max_refine'))
        if SC_type == SC_FF_CUB:
            SC = SC_freeform_Cubic(hi, ui, ki, max_refine)
        elif SC_type == SC_FF_BS:
            SC = SC_freeform_BSpline(hi, ui, ki, max_refine)
        elif SC_type == SC_FF_LIN:
            SC = SC_freeform_Linear(hi, ui, ki, max_refine)
    else:
        print('Unknown value of ''SC_type'': ', SC_type)
        exit(1)

    return SC


########################################################################
#              General (default) transformation functions              #
########################################################################

# Maximal value we allow to reach when optimizing unconstrained variable
TRANSFORM_MAX_VALUE = 1e50

def default_transformation(lbound, ubound, max_value = TRANSFORM_MAX_VALUE):
    """
        Transforms parameters on specified interval <ubound, lbound>
        and infinity values are replaces with 'max_value'.
        Transformation functions are chosen in such way that
        We consider 3 intervals:
        X1 = <lbound,   lbound+1>       Y1 = <-inf,     lbound+1>
        X2 = <lbound+1, ubound-1>       Y2 = <lbound+1, ubound-1>
        X3 = <ubound-1, ubound>         Y3 = <ubound-1, inf>

        We want function T(x), U(y) such that:
        T(X1) = Y1         U(Y1) = X1
        T(X2) = Y2         U(Y2) = X2
        T(X3) = Y3         U(Y3) = X3

             /-  lbound+1 + ln(x-lbound), x in X1
        T(x)=|-  x,                       x in X2
             \-  ubound-1 - ln(ubound-x), x in X3

             /-  lbound + e^(y-lbound-1), y in Y1
        U(y)=|-  y,                       y in Y2
             \-  ubound - e^(ubound-1-y), y in Y3

    """
    def transform_lb(x, lb):
        """ T(X1) """

        assert np.all(x > lb), \
            'Assertion x>lb failed.\nlb={}\nx={}'.format(lb, x)

        lb1 = lb + 1
        if np.isscalar(x):
            if x < lb1:
                x = lb1 + np.log(x - lb)
        else:
            for i in range(np.alen(x)):
                if x[i] < lb1:
                    x[i] = lb1 + np.log(x[i] - lb)

        return x

    def untransform_lb(y, lb):
        """ U(Y1) """

        lb1 = lb + 1;
        if np.isscalar(y):
            if y < lb1:
                y = lb + np.exp(y - lb1)
        else:
            for i in range(np.alen(y)):
                if y[i] < lb1:
                    y[i] = lb + np.exp(y[i] - lb1)

        return y

    def transform_rb(x, rb):
        """ T(X3) """

        assert np.all(x < rb), \
          'Assertion x<rb failed.\nrb={}\nx={}'.format(rb, x)

        rb1 = rb - 1
        if np.isscalar(x):
            if x > rb1:
                x = rb1 - np.log(rb - x)
        else:
            for i in range(np.alen(x)):
                if x[i] > rb1:
                    x[i] = rb1 - np.log(rb - x[i])

        return x

    def untransform_rb(y, rb):
        """ U(Y3) """
        rb1 = rb - 1;
        if np.isscalar(y):
            if y > rb1:
                y = rb - np.exp(rb1 - y)
        else:
            for i in range(np.alen(y)):
                if y[i] > rb1:
                    y[i] = rb - np.exp(rb1 - y[i])

        return y

    if (lbound == -np.inf) or (ubound == np.inf):
        lbound = max(lbound, -max_value)
        ubound = min(ubound, max_value)

        tr_fn = lambda x: transform_lb(transform_rb(x, ubound), lbound)
        un_fn = lambda y: untransform_lb(untransform_rb(y, ubound), lbound)
    else:
        I_half = (lbound + ubound)/2
        I_r    = (ubound-lbound)/np.pi

        tr_fn = lambda x: np.tan((x - I_half)/I_r)
        un_fn = lambda y: I_r * np.arctan(y) + I_half

    return (tr_fn, un_fn)
########################################################################
#                             Common class                             #
########################################################################

class SC_base():

    def __init__(self):
        self.h_init_max = 0.

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
                           h=None, rho=None, g=None, rtype='theta'):
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
            rtype    - return value type: 'theta', 'u', 'h', p

          Return values:
            theta    - saturation corresponding to pressure p /
                       relative saturation u / hydraulic head h
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
            h = -10.0* p /rho / g
            u = self.h2u(h)

        K     = self.u2Ku(u, Ks)

        if rtype == 'theta':
            xvalue = theta_r + (theta_s - theta_r) * u
        elif rtype == 'u':
            xvalue = u
        elif rtype in ('h', 'p'):
            if u is not None:
                h = self.u2h(u)

            if rtype == 'h':
                xvalue = h
            elif p is None:
                xvalue = - h * rho *g / 10 # p
            else:
                xvalue = p
        else:
            raise Exception('Unknown rtype: ', rtype)

        return (xvalue, K)

    def get_dyn_h_init(self, c_gammah, h_init_max):
        """
        Initial values of h~0 can cause troubles to the solver to start
        depending on parameters of the SC class. To ensure "smooth" start
        we compute a dynamically obtained 'h_init' value based on
        actual values. This may be important in the parameters
        optimization process.
        """
        self.h_init_max = h_init_max
        return self.h_init_max

    def add_transformations_fns(self, transform, untransform,
                                lbounds, ubounds):
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

    def canrefine_h(self):
        """
        indicate if this SC allows refinement of h
        """
        return False

    def refine(self, _):
        """
        refine the parameters and reinit. Return if success
        """
        return False

    def typeSC(self):
        """
        Indication of the compatible types of SC
        """
        raise NotImplementedError

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
        SC_base.__init__(self)
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
        self.h_init_max = min(c_gammah / self._gamma, h_init_max)
        return self.h_init_max

    def add_transformations_fns(self, transform, untransform,
                                lbounds, ubounds):
        """
        Transform/untransform methods for 'n' and 'gamma' parameters
        """
        max_value = TRANSFORM_MAX_VALUE

        transform['n']   = lambda n: max(np.log(n - 1.0), -max_value)
        untransform['n'] = lambda n_transf: 1+min(np.exp(n_transf), max_value)

        transform['gamma']   = lambda gamma: max(np.log(-gamma), -max_value)
        untransform['gamma'] = lambda gamma_transf: -min(np.exp(gamma_transf), max_value)

    def typeSC(self):
        """
        Indication of the compatible types of SC
        """
        return SC_vG

########################################################################
#                           Free-form  model                           #
########################################################################
class SC_freeform_base(SC_base):
    """
    A freefrom description of the saturation curve.
    In this, a discrete version of h is used: {h_i}, with -inf < h_i < 0.
    For these gridpoints, the value of effective saturation is passed:
        u_i = u(h_i)
    and also the value of the relative permeability:
        k_i = k(h_i)

    Depending on the interpolation used we can then determine u(h) and k(h)
    of the given points {(h_i, u_i)} and {(h_i, k_i)}.
    """
    def __init__(self, hi=None, ui=None, ki=None, refinemax=0,
                 compute_extra=False, issue_warning=False):
        SC_base.__init__(self)
        self.compute_extra = compute_extra
        self.extra = 0
        self.refinenr=0
        self.refinemax = refinemax
        self.oncheck_warning_only = issue_warning

        if hi is not None and ui is not None and ki is not None:
            SC_freeform_base.check_params(hi, ui, ki)
            self._set_values(np.array(hi, dtype=float), np.array(ui, dtype=float),
                              np.array(ki, dtype=float))
        else:
            self._hi = hi
            if np.iterable(hi):
                self._hi = np.array(self._hi)
            self._ui = ui
            self._ki = ki

        # NOTE: After base class initialization the subclass must run
        #       self.__interpolate_values() to get proper intializaton
        #       of the node values

    def _set_values(self, hi, ui, ki):
        #are internal parameters are such that all is monotone increasing
        # in terms of h
        self._hi = hi
        self._ui = ui
        self._ki = ki
        hmax = -min(1e-28, -hi[-1]*1e-3)
        if self.compute_extra:
            self.extra = 0
            if self.h_init_max > self._hi[-1] and self.h_init_max < hmax:
                self.extra = 1
        self._hnodes = np.empty(len(hi)+2+self.extra, float)
        self._uvals = np.empty(len(hi)+2+self.extra, float)
        self._kvals = np.empty(len(hi)+2+self.extra, float)
        #we use the log values  for head.
        self._hnodes[1:-1-self.extra] = -1*np.log(-hi[:])
        self._hnodes[0] = -1*np.log(-hi[0]*1e3)
        self._hnodes[-1] = -1*np.log(-hmax)
        self._hmax = -np.exp(-self._hnodes)
        #with log for head, saturation is linear
        self._uvals[1:-1-self.extra] = ui[:]
        self._uvals[0] = 0.
        self._uvals[-1] = 1.
        #with log for head, also rel perm should be log values
        self._kvals[1:-1-self.extra] = np.log(ki[:])
        self._kvals[0] = min(-50, self._kvals[1]-10 )
        self._kvals[-1] = 0.
        if self.compute_extra and (self.extra == 1):
            self._hnodes[-2] = -1*np.log(-self.h_init_max)
            self._uvals[-2]  = 1-1e-5
            self._kvals[-1]  = np.log(1-1e-5)

        # we now construct two cubic monotone interpolations: U(H) and K(H)
#        print ('test mono')
#        print ('hi', self._hi)
#        print ('ui', self._ui)
#        print ('ki', self._ki)
#
#        print ('hn', self._hnodes)
#        print ('uv', self._uvals)
#        print ('kv', self._kvals)

    def __interpolate_values(self):
        raise Exception ("Method must be subclassed.")

    @staticmethod
    def check_params(hi, ui, ki):
        if not np.iterable(hi) or not np.iterable(ui) or not np.iterable(ki):
            raise Exception ( 'Some parameters are not given: '\
                    'hi:%s, ui:%s, ki:%s' % (str(hi), str(ui), str(ki)))
        if not (len(hi)==len(ui)==len(ki)):
            raise Exception ('Parameters must have equal length')
        #hi must be monotone ascending > 0, and ui, ki monotone ascending
        ho = hi[0]
        uo = ui[0]
        ko = ki[0]
        for h,u,k in zip(hi[1:],ui[1:],ki[1:]):
            if not h>ho or not h<0.:
                raise Exception('Hydraulic head h must be negative and a '
                                'monotone ascending array, instead %s' % str(hi))
            if not u>uo or not uo>0:
                raise Exception('Effective saturation Se must be positive and a '
                                'monotone ascending array in terms of h, instead %s' % str(ui))
            if not k>ko or not ko>0:
                raise Exception('Relative permeability k must be positive and a '
                                'monotone ascending array in terms of h, instead %s' % str(ki))

        #The edges of ui and ki are fixed and will not be optimized
        if not (ui[0] >0) or not (ui[-1] < 1):
            if self.oncheck_warning_only:
                print ("WARNING: Effective saturation Se starts at 0, ends "
                       "at 1, only pass values between these extremes!")
            else:
                raise Exception("Effective saturation Se starts at 0, ends at "
                                "1, only pass values between these extremes!")
        if not (ki[0] > 0) or not (ki[-1] < 1):
            raise Exception('Relative permeability k starts at 0, ends at 1, '
                            'only pass values between these extremes!')

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
        SC_freeform_base.check_params(self._hi, ui, ki)
        self._set_values(self._hi, ui, ki)

    def get_parameters(self):
        """
        return the current parameters
        """
        return self._ui, self._ki

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
        wherezero = (h >= 0)
        wheresmall = (h < -np.exp(-self._hnodes[0]))
        tmp1 = self.logh2u(-np.log(-h) )
        if not u is None:
            u[:] = tmp1[:]
        else:
            u    = tmp1
        u[wherezero] = 1.
        #less than smallest in our approx, set u to zero
        u[wheresmall] = 0.
        return u

    def u2h(self, u, h = None):
        """
        Return the value of h for the given effective saturation u and
        parameters passed
        """
        internal_u = u[:]
        internal_u[ u<self._uvals[ 0] ] = self._uvals[ 0]
        internal_u[ u>self._uvals[-1] ] = self._uvals[-1]
        tmph = self.logh2u.root(internal_u)
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

    def add_transformations_fns(self, transform, untransform,
                                lbounds, ubounds):
        """
        Transform/untransform methods for 'ki' and 'ui' parameters
        """
        # ui, ki are on <0, 1>
        #(transform['ki'], untransform['ki']) = \
        #  default_transformation(lbounds['ki'], ubounds['ki'])

        #(transform['ui'], untransform['ui']) = \
        #  default_transformation(lbounds['ui'], ubounds['ui'])

        transform['ui']   = lambda ui: ui
        untransform['ui'] = lambda ui_transf: ui_transf

        transform['ki']   = lambda ki: np.log(ki)
        untransform['ki'] = lambda ki_transf: np.exp(ki_transf)


class SC_freeform_Cubic(SC_freeform_base):
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
    def __init__(self, hi=None, ui=None, ki=None, refinemax=0):
        SC_freeform_base.__init__(self, hi, ui, ki, refinemax,
                                  compute_extra=True, issue_warning=False)
        self.__interpolate_values()


    def __interpolate_values(self):
        self.logh2u = MonoCubicInterp(self._hnodes, self._uvals)
        self.logh2logk = MonoCubicInterp(self._hnodes, self._kvals)

    def get_dyn_h_init(self, c_gammah, h_init_max):
        """
        Initial values of h~0 can cause troubles to the solver to start
        depending on parameters of the SC class. To ensure "smooth" start
        we compute a dynamically obtained 'h_init' value based on
        actual values. This may be important in the parameters
        optimization process.
        """
        recreate = False
        if self._hi is not None and self.h_init_max != h_init_max:
            #we need to regenerate the bsplace!
            recreate = True
        self.h_init_max = h_init_max
        if recreate and self._hi is not None and self._ui is not None \
                and self._ki is not None:
            self._set_values(self._hi, self._ui, self._ki)

        return self.h_init_max

    def typeSC(self):
        """
        Indication of the compatible types of SC
        """
        return SC_FF_CUB

    def canrefine_h(self):
        """
        indicate if this SC allows refinement of h
        """
        return True

    def refinable_h(self):
        """
        Obtain the current h values which are allowed to be refined
        """
        return -np.exp(-self._hnodes[1:-1])

    def refine(self, prev_measurements):
        """
        refine the parameters and reinit. Return if success
        """
        if self.refinenr >= self.refinemax:
            return False
        #we refine the parameter space. Here we add a h point, and hence
        #add 2 new parameters: the ki and ui at that new hi.
        # obtain (x, h) computed data and (x, u)
        comph  = prev_measurements.get_computed_value('h')[1]
        compu  = prev_measurements.get_computed_value('u')[1]
        compui = self._ui
        #compki = self._ki
        if self.refinenr == 0:
            #first refine. refine strategy must start from 2 points!
            self._internalrefine = 1
            if not len(self._hi) == 2:
                raise Exception("Refinement must start from minimum amount "
                                "of points: 2 given head values")

        hn = self._hnodes[1:-1-self.extra]
        minh = np.min(comph.flatten())

        while True:
            if self._internalrefine > 32:
                #stop computation:
                print ("More than maximum number of internal refinements!")
                return False
            curexp = 1+int(np.trunc(np.log2(self._internalrefine)))
            curintervalsize = 1/2**(curexp)
            print ('intref', self._internalrefine)
            possiblenewhnodes = [hn[0] + (2*i+1)*curintervalsize * (hn[-1]-hn[0])
                                                for i in range(2**(curexp-1))]
            print ('newhn 0', possiblenewhnodes, 'h=', -np.exp(-np.array(possiblenewhnodes)))
            possiblenewhnodes = [i for i in possiblenewhnodes if not i in hn]
            print ('newhn 1', possiblenewhnodes)
            #we select only h > minh
            possiblenewhnodes = np.array([i for i in possiblenewhnodes
                                            if -np.exp(-i) > minh*2])
            print ('newhn E', possiblenewhnodes, 'h=', -np.exp(-np.array(possiblenewhnodes)))
            if len(possiblenewhnodes) == 0:
                self._internalrefine += 1
                continue
            else:
                possiblenewh = -np.exp(-possiblenewhnodes)
                #start with the one in interval with largest u change
                uchange = 0;
                prevui = self._uvals[1]
                nexth = possiblenewh[0]
                nexthnode = possiblenewhnodes[0]
                for h, hnode in zip(possiblenewh, possiblenewhnodes):
                    for hi, ui in zip(self._hi, self._ui):
                        if hi > h:
                            diffu = ui-prevui
                            if diffu > uchange:
                                nexth = h
                                nexthnode = hnode
                                uchange = diffu
                        prevui = ui
                break;

        #we have a new h to add, we insert it, and update parameters
        nh = self.h2u(np.array([nexth]))
        nk = self.h2Kh(np.array([nexth]), 1)
        ind = np.searchsorted(self._hnodes, nexthnode)
        self._hnodes = np.insert(self._hnodes, ind, nexthnode)
        self._hi     = np.insert(self._hi, ind-1, nexth)
        self._ui     = np.insert(self._ui, ind-1, nh)
        self._ki     = np.insert(self._ki, ind-1, nk)
        self._uvals  = np.insert(self._uvals, ind, self._ui[ind-1])
        self._kvals  = np.insert(self._kvals, ind, np.log(self._ki[ind-1]))

        print ("*** refined SC curve  - new par ***")
        print (" h  ", self._hi)
        print (" u  ", self._ui)
        print (" k  ", self._ki)
        print (" hn ", self._hnodes)
        print (" uv ", self._uvals)
        print (" kv ", self._kvals)

        self.__interpolate_values()

        self.refinenr += 1
        self._internalrefine += 1
        return True

########################################################################
#                  Free-form  model local BSpline based                #
########################################################################

class SC_freeform_BSpline(SC_freeform_base):
    """
    A freefrom description of the saturation curve.
    In this, a discrete version of h is used: {h_i}, with -inf < h_i < 0.
    For these gridpoints, the value of effective saturation is passed:
        u_i = u(h_i)
    and also the value of the relative permeability:
        k_i = k(h_i)

    We then determine u(h) and k(h) as a localized quadratic B-spline
    of the given points {(h_i, u_i)} and {(h_i, k_i)}
    From this cubic spline the derivatives can be extracted
    """
    def __init__(self, hi=None, ui=None, ki=None, refinemax=0):
        SC_freeform_base.__init__(self, hi, ui, ki, refinemax,
                                  compute_extra=False, issue_warning=True)
        self.__interpolate_values()


    def __interpolate_values(self):
        self.logh2u = QuadraticBspline(self._hnodes, self._uvals)
        self.logh2logk = QuadraticBspline(self._hnodes, self._kvals)

    def typeSC(self):
        """
        Indication of the compatible types of SC
        """
        return SC_FF_BS

class SC_freeform_Linear(SC_freeform_base):
    def __init__(self, hi=None, ui=None, ki=None, refinemax=0):
        SC_freeform_base.__init__(self, hi, ui, ki, refinemax,
                                  compute_extra=True, issue_warning=False)
        self.__interpolate_values()

    def __interpolate_values(self):
        self.logh2u    = PiecewiseLinear(self._hnodes, self._uvals)
        self.logh2logk = PiecewiseLinear(self._hnodes, self._kvals)

    def typeSC(self):
        return SC_FF_LIN

    def canrefine_h(self):
        """
        indicate if this SC allows refinement of h
        """
        return True

    def refine(self, prev_measurements):

        if self.refinenr >= self.refinemax:
            return False

        # Probably not needed for piecewise linear interpolation
        # #first refine must start from 2 points
        # if (self.refinenr == 0) and (if not len(self._hi) == 2):
        #     raise Exception("Refinement must start from minimum amount "
        #                     "of points: 2 given head values")

        # Copy & paste of the CubicInterpolator.refine()
        hn = self._hnodes[1:-1-self.extra]
        minh = np.min(comph.flatten())

        while True:
            if self._internalrefine > 32:
                #stop computation:
                print ("More than maximum number of internal refinements!")
                return False

            curexp = 1+int(np.trunc(np.log2(self._internalrefine)))
            curintervalsize = 1/2**(curexp)
            print ('intref', self._internalrefine)
            possiblenewhnodes = [hn[0] + (2*i+1)*curintervalsize * (hn[-1]-hn[0])
                                                for i in range(2**(curexp-1))]
            print ('newhn 0', possiblenewhnodes)
            possiblenewhnodes = [i for i in possiblenewhnodes if not i in hn]
            print ('newhn 1', possiblenewhnodes)
            #we select only h > minh
            possiblenewhnodes = np.array([i for i in possiblenewhnodes
                                            if -np.exp(-i) > minh])
            print ('newhn E', possiblenewhnodes)
            if len(possiblenewhnodes) == 0:
                self._internalrefine += 1
                continue
            else:
                possiblenewh = -np.exp(-possiblenewhnodes)
                #start with the one in interval with largest u change
                uchange = 0;
                prevui = self._uvals[1]
                nexth = possiblenewh[0]
                nexthnode = possiblenewhnodes[0]
                for h, hnode in zip(possiblenewh, possiblenewhnodes):
                    for hi, ui in zip(self._hi, self._ui):
                        if hi > h:
                            diffu = ui-prevui
                            if diffu > uchange:
                                nexth = h
                                nexthnode = hnode
                                uchange = diffu
                        prevui = ui
                break;

        #we have a new h to add, we insert it, and update parameters
        nh = self.h2u(np.array([nexth]))
        nk = self.h2Kh(np.array([nexth]), 1)
        ind = np.searchsorted(self._hnodes, nexthnode)
        self._hnodes = np.insert(self._hnodes, ind, nexthnode)
        self._hi     = np.insert(self._hi, ind-1, nexth)
        self._ui     = np.insert(self._ui, ind-1, nh)
        self._ki     = np.insert(self._ki, ind-1, nk)
        self._uvals  = np.insert(self._uvals, ind, self._ui[ind-1])
        self._kvals  = np.insert(self._kvals, ind, np.log(self._ki[ind-1]))

        print ("*** refined SC curve  - new par ***")
        print (" h  ", self._hi)
        print (" u  ", self._ui)
        print (" k  ", self._ki)
        print (" hn ", self._hnodes)
        print (" uv ", self._uvals)
        print (" kv ", self._kvals)

        self.__interpolate_values()

        self.refinenr += 1
        self._internalrefine += 1
        return True

if __name__ == "__main__":
    #test of SC curves
    #BSPLINE on hi, ui, ki
    ks = 0.00011846618
    hi = [-403.43, -200, -1]
    ui = [0.73394395, 0.8625, 0.9923954]
    ki = np.array([0.00239603,  0.043, 0.08027015])
    #BSPLINE on hi, ui, ki
    hiext = [-403.43*1e3, -403.43, -200, -1, -1*1e-3]
    uiext = [0          , 0.73394395,  0.8625, 0.9923954, 1]
    kiext = [0.00239603*np.exp(-10), 0.00239603, 0.043/ks, 0.08027015, 1]
    BS = SC_freeform_BSpline(hi, ui, ki)
    FF = SC_freeform_Cubic(hi, ui, ki)
    LI = SC_freeform_Linear(hi, ui, ki)

    haxis = np.linspace(-500., -0.0001, 10000)
    uax = BS.h2u(haxis)
    kax = BS.h2Kh(haxis, ks)
    duax = BS.dudh(haxis)
    uaxf = FF.h2u(haxis)
    kaxf = FF.h2Kh(haxis, ks)
    duaxf = FF.dudh(haxis)
    uaxl  = LI.h2u(haxis)
    kaxl  = LI.h2Kh(haxis,ks)
    duaxl = LI.dudh(haxis)

    import pylab
    pylab.figure(1)
    pylab.xlim(xmin=hi[0])
    pylab.title("effective saturation")
    pylab.plot(hiext, uiext, 'k-')
    pylab.plot(haxis, uax, 'b-', label="BS")
    #pylab.plot(haxis, duax, 'r-')
    pylab.plot(haxis, uaxf, 'g-', label="FF")
    #pylab.plot(haxis, duaxf, 'c-')
    pylab.plot(haxis, uaxl, 'r-', label="LF")
    pylab.legend()
    pylab.show()

    pylab.figure(2)
    pylab.title("relative permeability")
    pylab.plot(hi, ki, 'k-')
    pylab.plot(haxis, kax/ks, 'b-', label="BS")
    pylab.plot(haxis, kaxf/ks, 'g-', label="FF")
    pylab.plot(haxis, kaxl/ks, 'r-', label="LF")
    pylab.ylim(ymax=1.1)
    pylab.legend()
    pylab.show()
