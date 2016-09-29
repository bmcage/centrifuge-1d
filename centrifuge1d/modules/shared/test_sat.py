# -*- coding: utf-8 -*-
"""
Created on Thu Sep 29 12:02:19 2016

@author: benny
"""
from __future__ import print_function, division

import numpy as np
from scipy.optimize import bisect, minimize
NEW_LS = True
try:
    from scipy.optimize import least_squares
except:
    NEW_LS = False
    from scipy.optimize import leastsq as least_squares

import pylab

from interpolate import (MonoCubicInterp, QuadraticBspline, PiecewiseLinear,
                         PiecewiseLinearMonotoneAsc)
from saturation_curve import SC_Durner_freeform

__NR = 400
P_DEFAULT = np.linspace(-1, 9, __NR)
P_DEFAULT = np.power(10* np.ones(__NR), P_DEFAULT)
P_DEFAULT[0] = 0

p = P_DEFAULT[1:]
rho = 1.0
g = 981.
h = -10.0* p /rho / g
test_logh =  -1*np.log(-h)

"""
20p NP fines:
Optimal parameters found:
  gamma2   : [-0.01198936]
  gamma1   : [-0.6459756]
  ki       : [  2.55592585e-177   8.37197116e-084   2.45017122e-077   2.41964065e-004
   1.89611978e-003   4.73289371e-001]
  Ks [cm/s]: [  7.91674655e-07]
  w1       : [ 0.17095876]
  hiadd    : [-69.04066976]
  n1       : [ 47.21602695]
  n2       : [ 2.01094139]

Used parameters:
  Ks [cm/s]: [  7.91674655e-07]
  hi       : [-800.         -184.69972153  -69.04066976  -84.88393601  -25.55339634
   -1.        ]

(h,S) values of experiment:
0.01	1
1.5720209518	0.8853685552
14.013437767	0.8255608449
55.6122655858	0.6660736173
124.9067220385	0.461730607
498.7431928887	0.1577080795
885.4760794438	0.0630125381
1198.8125015337	0.0131727795

We determine ki values at the same h values
"""
TEST = '20pNPF'

if TEST == '20pNPF':
    gamma2   = -0.01198936
    gamma1   = -0.6459756
    Ks       = 7.91674655e-07
    w1       = 0.17095876
    n1       = 47.21602695
    n2       = 2.01094139
    #hiorig =  [-800.,  -184.69972153,  -69.04066976,  -84.88393601,  -25.55339634, -1.]
    hiorig = [-800., -184.69972153, -84.88393601, -25.55339634, -1.]
    #kiorig = [  2.55592585e-177, 8.37197116e-084, 2.45017122e-077,   2.41964065e-004,
    #            1.89611978e-003,   4.73289371e-001]
    kiorig = [  4.28434076e-177,   1.54026912e-083,   7.29970690e-005,
                2.33028423e-003,   7.03968739e-001]

    #hinew = [-1198.8, -885.5, -498.7, -124.9, -55.6, -14.0, -1.6, -0.01]
    hinew = [-1198.8, -885.5, -498.7, -124.9, -110, -90, -70, -55.6, -14.0, -1.6, -0.01]

else:
    raise ValueError, 'Unknown test %s' % TEST

hiorig = np.asarray(hiorig)
kiorig = np.asarray(kiorig)
hinew = np.asarray(hinew)

hnodesorig = np.empty(len(hiorig)+2, float)
kvalsorig = np.empty(len(hiorig)+2, float)
#we use the log values  for head.

hmax = -min(1e-28, -hiorig[-1]*1e-3)
hnodesorig[1:-1] = -1*np.log(-hiorig[:])
hnodesorig[0] = -1*np.log(-hiorig[0]*1e3)
hnodesorig[-1] = -1*np.log(-hmax)
#with log for head, also rel perm should be log values
kvalsorig[1:-1] = np.log(kiorig[:])
kvalsorig[0] = min(-50, kvalsorig[1]-10 )
kvalsorig[-1] = 0.

print ('Given hi:', hiorig, 'with ki', kiorig, hnodesorig.shape)
print ('Given hnodes:', hnodesorig, 'with kvals', kvalsorig)
logh2logk = MonoCubicInterp(hnodesorig, kvalsorig)

hnodesnew = np.empty(len(hinew)+2, float)
kvalsnew = np.empty(len(hinew)+2, float)
hnodesnew[1:-1] = -1*np.log(-hinew[:])
hnodesnew[0] = -1*np.log(-hiorig[0]*1e3)
hnodesnew[-1] = -1*np.log(-hmax)
kvalsnew = logh2logk(hnodesnew)

def _mink_fun_res(x):
    """ internally used function to minimize in order to find new kval
        after refine. This function computes the residuals"""
    new_kvals = np.empty(len(x)+2, float)
    new_kvals[1:-1] = x[:]
    new_kvals[0] = kvalsorig[0]
    new_kvals[-1] = kvalsorig[-1]

    # kvals must be monotone increasing, check this:
    ko = x[0]
    error =0.
    for k in x[1:]:
        if (not k>ko):
            print ('ERROR in checking k>ko', k, ko)
            error += (ko*1.001 - k)*1e5
        ko = k

    #The edges of ki are fixed and will not be optimized
    if not (x[0] > 0):
        error += (-x[0]+1e-3)*1e5
    if not (x[-1] < 1):
        error +=  (x[-1]-1+1e-3)*1e5

    new_logh2logk = MonoCubicInterp(hnodesnew, new_kvals)
    return np.abs(new_logh2logk(test_logh)
            -  logh2logk(test_logh)) + error/len(x)


res_minkval = least_squares(_mink_fun_res, kvalsnew[1:-1])
if NEW_LS:
    if res_minkval.success:
        print('Found optimial start after a refine in', res_minkval.x,
              'Startval was', kvalsnew[1:-1])
    else:
        logmessage = 'Refine FAILED, msg=' + res_minkval.message
        print(logmessage)
        raise Exception(res_minkval.message)
    kvalsnew[1:-1] = res_minkval.x
else:
    print ('Full output refine least_squares:', res_minkval)
    kvalsnew[1:-1] = res_minkval[0]

#we test the result
ko = kvalsnew[0]
_REFINEERROR = False
for k in kvalsnew[1:]:
    if (not k>ko):
        print ('ERROR in checking k>ko', k, ko)
        _REFINEERROR = True
    ko = k
if _REFINEERROR:
    #resort to avoid problems ....
    print ('orig kval', kvalsnew)
    kvalsnew[:] = sorted(kvalsnew)
    print ('changed kval', kvalsnew)

#adapt ki based on determined kvals:
kinew = np.exp(kvalsnew[1:-1] )

print ('Found solution')
print ('hi new nodes =', hinew)
print ('start ki values =', kinew)


SWRCorig = SC_Durner_freeform(n1=n1, gamma1=gamma1, n2=n2, gamma2=gamma2,
                          w1=w1, hi=hiorig, ki=kiorig)
SWRC = SC_Durner_freeform(n1=n1, gamma1=gamma1, n2=n2, gamma2=gamma2,
                          w1=w1, hi=hinew, ki=kinew)


haxis = np.linspace(-500., -0.0001, 10000)
uax = SWRC.h2u(haxis)
kax = SWRC.h2Kh(haxis, Ks)
uaxorig = SWRCorig.h2u(haxis)
kaxorig = SWRCorig.h2Kh(haxis, Ks)

pylab.figure(1)
#pylab.ylim(ymax=-hinew[0])
pylab.title("effective saturation")
pylab.ylabel('Negative Head [cm]')
pylab.yscale('log')
pylab.xlabel('Saturation')
pylab.plot(uaxorig, -haxis, 'k-', label="Orig")
pylab.plot(uax, -haxis, 'b-', label="New")
pylab.legend()


pylab.figure(2)
pylab.title("Hydraulic Conductivity")
pylab.ylabel('Hydraulic Conductivity [cm/s]')
pylab.yscale('log')
pylab.xlabel('Saturation')
pylab.ylim(ymin=1e-20)
pylab.plot(uaxorig, kaxorig, 'k-', label="Orig")
pylab.plot(uax, kax, 'b-', label="New")
pylab.legend()

pylab.figure(3)
pylab.title("Hydraulic Conductivity")
pylab.ylabel('Hydraulic Conductivity [cm/s]')
pylab.xlabel('Negative Head [cm]')
pylab.yscale('log')
pylab.xscale('log')
pylab.ylim(ymin=1e-40)
pylab.plot(-haxis, kaxorig, 'k-', label="Orig")
pylab.plot(-haxis, kax, 'b-', label="New")
pylab.plot(-hinew, kinew*Ks, 'kx', label="Grid")
pylab.legend()

pylab.show()
