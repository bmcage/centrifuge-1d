from __future__ import division, print_function

"""
  This modules contains consolidation object(s) for single flow
"""
import numpy as np
from interpolate import (MonoCubicInterp, QuadraticBspline, PiecewiseLinear,
                         PiecewiseLinearMonotoneAsc)
from saturation_curve import (TRANSFORM_MAX_VALUE,
                              default_transformation)


__NR = 400
E_DEFAULT = np.linspace(0.01, 9, __NR)

CON_SLURRY     = 1
CON_GOMPERTZ   = 2
CON_FREEFORM   = 3

def get_dict_value(dictionary, keys):
    if type(keys) is str:                 # single item
        return dictionary[keys]
    else:                                 # list of keys
        return [dictionary[key] for key in keys]

def create_CON(data):
    """
    Returns a new CON object. Data is expected to be either a dict with values
    or an instance of Configuration class.
    Values present:
        {'con_type': ("Type of constitutive law to use for consolidation "
                        "sigma(e) and K(e). slurry=1, preconsolidated=2, freeform=3"),
         'con_max_refine': ("For con_type freeform with refinement, how many "
                         "times are we allowed to refine the grid?"),
         'a' : ("A parameter for the constitutive consolidation laws"),
         'b' : ("B parameter for the constitutive consolidation laws"),
         'c' : ("C parameter for the constitutive consolidation laws"),
         'd' : ("D parameter for the constitutive consolidation laws"),
         'cc': ("Cc parameter, the compression index"),
         'ei' : ("For freeform consolidation grid points in void radio e"),
         'si' : ("For freeform consolidation effective stress sigma values "
                 "with the grid points in void radio e"),
         'ki' : ("For freeform consolidation saturated conductivity values "
                 "with the grid points in void radio e"),
         'eiadd' : ("For freeform consolidation where to add the next "
                     "refinement point in grid points in void radio e"),
         'e0' : ("Initial void ratio")
        },
    """
    if type(data) is dict:                # dict
        get_value = lambda keys: get_dict_value(data, keys)
    else:                                 # Configuration class instance
        get_value = data.get_value

    CON_type = get_value('con_type')
    if CON_type == CON_SLURRY:
        (e0, A, B, C, D) = get_value(('e0', 'a', 'b', 'c', 'd'))
        CON = CON_Slurry(e0, A, B, C, D)
    elif CON_type == CON_GOMPERTZ:
        (e0, A, B, Cc, C, D) = get_value(('e0', 'a', 'b', 'cc', 'c', 'd'))
        CON = CON_Gompertz(e0, Cc, A, B, C, D)
    else:
        print('Unknown value (or not yet implemented) of ''CON_type'': ', CON_type)
        exit(1)

    return CON


########################################################################
#                             Common class                             #
########################################################################

class CON_base():
    """
    Base class for conductivity relationships
    Note:
    sigmaprime is effective stress in unit g/(s^2 cm) = 0.1 Pa. Typical
        consolidation relationships are represented in kPa, values from 1 to 100
    K is hydraulic conductivity in cm/s
    """

    def __init__(self):
        pass

    def add_transformations_fns(self, transform, untransform,
                                lbounds, ubounds):
        """
        When perfoming optimization problem over the CON parameters,
        we can work directly with the parameters values or
        some transformation of them (e.g. to ensure we keep the
         searching interval bounded).
        This parameters transformation is performed ONLY to determine
        next value in the inverse optimization process, i.e.methods
        using these parameters always obtain untransformed values.
        """
        pass

    def canrefine_e(self):
        """
        indicate if this CON allows refinement of e
        """
        return False

    def refine(self, _1, _2, _3, _4, _5):
        """
        refine the parameters and reinit. Return if success
        """
        return False

    def typeCON(self):
        """
        Indication of the compatible types of SC
        """
        raise NotImplementedError

    def effstress_curve(self, e=None):
        """
          Determine the effective stress curve.

          Parameters:
            e        - void ratio

          Return values:
            e                 - void ratio
            effective stress  - effective stress sigma^' for the e values
        """
        if (e is None):
            e = E_DEFAULT

        return (e, self.e2sigmaprime(e))


    def conductivity_curve(self, e=None):
        """
          Determine the conductivity curve.

          Parameters:

          Return values:
            e        - void ratio
            K        - saturated hydraulic conductivity
        """

        if (e is None):
            e = E_DEFAULT

        K     = self.e2Ks(e)

        return (e, K)

########################################################################
#                       Slurry based model                             #
########################################################################

class CON_Slurry(CON_base):
    """ Slurry model.
    We start consolidation with void ratio e0 being the fluid limit
    This should be good model for high e, so e from eg 3 to 7

    We have
        effective stress
         sigprime(e) = A ((e0-e)/(1+e0))^B
            A,B constants, with B close to 2
        hydraulic conductivity
         K(e)  =  (1+e) (C+De)

    Typical values: A := 2e6; B := 1.8; e0 := 7.; C := 0.1e-6; D := 0.13e-5
    """
    def __init__(self, e0=None, A=None, B=None, C=None, D=None):
        CON_base.__init__(self)
        self._e0 = e0
        self._A = A
        self._B = B
        self._C = C
        self._D = D

    def set_parameters(self, params):
        if 'e0' in params:
            self._e0 = params['e0']

        if 'a' in params:
            self._A = params['a']
        if 'b' in params:
            self._B = params['b']
        if 'c' in params:
            self._C = params['c']
        if 'd' in params:
            self._D = params['d']

    def e2Ks(self, e, Ks = None):
        """
        For a given void ratio e, what is the hydraulic conductivity
        K(e)  =  (1+e) (C+De)
        """
        if not Ks is None:
            Ks[:] = (1+e)*(self._C + self._D * e)
        else:
            Ks = (1+e)*(self._C + self._D * e)
        return Ks

    def e2sigmaprime(self, e, sigp = None):
        """
        For a given void ratio, what is the sigma prime (effective stress)
        sigprime(e) = A ((e0-e)/(1+e0))^B
        """
        einternal = np.empty(len(e), float)
        einternal[:] = e[:]
        #numerically, we can have values of e>e0, we correct for those
        einternal[self._e0 < e] = self._e0
        if not sigp is None:
            sigp[:] = self._A  * np.power((self._e0 - einternal) / (1 + self._e0), self._B)
        else:
            sigp = self._A  * np.power((self._e0 - einternal) / (1 + self._e0), self._B)
        return sigp

    def sigmaprime2e(self, sigp, e = None):
        """
        For a given value of of sigma prime, what is the void ratio

        """
        if not e is None:
            e[:] = self._e0 - (1+self._e0) * np.power(sigp/self._A, 1/self._B)
        else:
            e = self._e0 - (1+self._e0) * np.power(sigp/self._A, 1/self._B)
        return e

    def dKsde(self, e, dKsde = None):
        """
        Return value of d Ks/ de in e
        """
        if not dKsde is None:
            dKsde[:] = self._D * (1+2*e) + self._C
        else:
            dKsde = self._D * (1+2*e) + self._C
        return dKsde

    def dKo1pede(self, e, dKo1pede=None):
        """
        Return value of d (K/(1+e)) / de in e
        """
        if not dKo1pede is None:
            dKo1pede[:] = self._D
        else:
            dKo1pede = 0*e+self._D
        return dKo1pede

    def dsigpde(self, e, dsigpde = None, zeroval=0.):
        """
        return value of d sigma' / de at e
        Note: for CON_Slurry, values at e0 go to 0, so we avoid derivatives
              there, enforcing e<=0.999 e0. This is a needed regularization
        If zeroval != 0, all values > zeroval are given the value zeroval.
          For example: zeroval = -1e-10 returns as highest value -1e10!
        """
        einternal = np.empty(len(e), float)
        einternal[:] = e[:]
        #numerically, we can have values of e>e0, we correct for those
        einternal[self._e0*0.999 < einternal] = self._e0*.999 #2*self._e0-e[self._e0<e]
        if not dsigpde is None:
            dsigpde[:] = -self._B * self._A  \
                * np.power((self._e0 - einternal), self._B-1) \
                / np.power((1 + self._e0), self._B)
        else:
            dsigpde = -self._B * self._A  \
                * np.power((self._e0 - einternal), self._B-1) \
                / np.power((1 + self._e0), self._B)
        if zeroval != 0:
            dsigpde[dsigpde>zeroval] = zeroval
        return dsigpde

    def typeCON(self):
        """
        Indication of the compatible types of CON
        """
        return CON_SLURRY


########################################################################
#                       Reconsolidation based on Gompertz function     #
########################################################################

class CON_Gompertz(CON_base):
    """ Reconsolidation Gompertz model.
    We start consolidation with uniform void ratio e0
    For loading (consolistion) we use Gompertz, e0, Cc
    This should be good model for high e, so e from eg 3 to 7

    We have
        effective stress
         sigprime(e) = a+c*exp(-exp(b*(log10(effstress)-m)))
            We use c=A, m=B, a=e0-c, b=exp(1)*Cc/c;
            A,B constants, e0 initial void ratio, Cc compression index
        hydraulic conductivity
         K(e)  =  (1+e) (C+De)

    Typical values: Cc := 1.2; c=A := 3.1; e0 := 4.; m=B := log10(20000) = 4.301;
                    C := 0.1e-6; D := 0.13e-5
    """
    def __init__(self, e0 = None, Cc=None, A=None, B=None, C=None, D=None):
        CON_base.__init__(self)
        self._e0 = e0
        self._Cc = Cc
        self._c = A
        self._m = B
        self._C = C
        self._D = D

    def set_parameters(self, params):
        if 'e0' in params:
            self._e0 = params['e0']

        if 'cc' in params:
            self._Cc = params['cc']
        if 'a' in params:
            self._c = params['a']
        if 'b' in params:
            self._m = params['b']
        if 'c' in params:
            self._C = params['c']
        if 'd' in params:
            self._D = params['d']

    def e2Ks(self, e, Ks = None):
        """
        For a given void ratio e, what is the hydraulic conductivity
        K(e)  =  (1+e) (C+De)
        """
        if not Ks is None:
            Ks[:] = (1+e)*(self._C + self._D * e)
        else:
            Ks = (1+e)*(self._C + self._D * e)
        return Ks

    def e2sigmaprime(self, e, sigp = None):
        """
        For a given void ratio, what is the sigma prime (effective stress)
        sigprime(e) =
        """
        einternal = np.empty(len(e), float)
        einternal[:] = e[:]
        #numerically, we can have values of e>e0, we correct for those
        if len(e) == 1:
            einternal[0] = self._e0
            if einternal[0]> self._e0: einternal[0]=self._e0
        else:
            einternal[self._e0 < e] = self._e0
        #note, einternal > a is required !!
        a = self._e0 - self._c
        b = np.exp(1)*self._Cc / self._c
        if not sigp is None:
            sigp[:] = np.exp((np.exp(1)*self._Cc*self._m  \
                              + np.log(np.log(self._c/(einternal-a)))*self._c) \
                             *np.log(10)/(self._Cc*np.exp(1)))
        else:
            sigp = np.exp((np.exp(1)*self._Cc*self._m  \
                              + np.log(np.log(self._c/(einternal-a)))*self._c) \
                             *np.log(10)/(self._Cc*np.exp(1)))
        #print ('test', self._Cc, np.log(self._c/(einternal-a)))
        #raw_input(' test ' + str(sigp) )
        return sigp

    def sigmaprime2e(self, sigp, e = None):
        """
        For a given value of of sigma prime, what is the void ratio

        """
        a = self._e0 - self._c
        b = np.exp(1)*self._Cc / self._c
        if not e is None:
            e[:] = a+self._c*np.exp(-np.exp(b*(np.log(sigp)/np.log(10)-self._m)))
        else:
            e = a+self._c*np.exp(-np.exp(b*(np.log(sigp)/np.log(10)-self._m)))
        return e

    def dKsde(self, e, dKsde = None):
        """
        Return value of d Ks/ de in e
        """
        if not dKsde is None:
            dKsde[:] = self._D * (1+2*e) + self._C
        else:
            dKsde = self._D * (1+2*e) + self._C
        return dKsde

    def dKo1pede(self, e, dKo1pede=None):
        """
        Return value of d (K/(1+e)) / de in e
        """
        if not dKo1pede is None:
            dKo1pede[:] = self._D
        else:
            dKo1pede = 0*e+self._D
        return dKo1pede

    def dsigpde(self, e, dsigpde = None, zeroval=0.):
        """
        return value of d sigma' / de at e

        If zeroval != 0, all values > zeroval are given the value zeroval.
          For example: zeroval = -1e-10 returns as highest value -1e10!
        """
        einternal = np.empty(len(e), float)
        einternal[:] = e[:]
        #numerically, we can have values of e>e0, we correct for those
        einternal[self._e0*0.999 < einternal] = self._e0*.999 #2*self._e0-e[self._e0<e]
        sigp = self.e2sigmaprime(e, None)
        a = self._e0 - self._c
        if not dsigpde is None:
            dsigpde[:] = - self._c * np.log(10)*sigp \
                / ((einternal-a)*np.log(self._c/(einternal-a))
                                *self._Cc*np.exp(1)
                  )
        else:
            dsigpde = - self._c * np.log(10)*sigp \
                / ((einternal-a)*np.log(self._c/(einternal-a))
                                *self._Cc*np.exp(1)
                  )

        #print ('dsipde', dsigpde)
        #raw_input('test')
        #if zeroval != 0:
        #    dsigpde[dsigpde>-1e-8] = -1e-8
        return dsigpde

    def typeCON(self):
        """
        Indication of the compatible types of CON
        """
        return CON_GOMPERTZ