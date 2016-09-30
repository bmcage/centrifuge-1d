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
from saturation_curve import (SC_Durner, SC_Durner_freeform, SC_vanGenuchten)

STYLES = ['k-', 'b-', 'r-', 'c-', 'g-']

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

REFINE = False
#TESTS = ['20pNPF_DurnerFF_try4', '20% NP fines', '10% fines']
TESTS = ['galaa-20pNPfines-dr_DurnerFF_try5-002', '20% NP fines',
         'galaa-20pfines-dr_DurnerFF_try6-002',
         '30% NP fines', '20% fines', '10% NP fines',
         'galaa-10pfines-dr_DurnerFF_try6-006',
         '10% fines', 'galaa-15pfines-dr_DurnerFF_try7-002',
         '20pNPF_DurnerFF_try4', 'galaa-30pNPfines-dr_DurnerFF_try6-002',
         '15% fines', 'galaa-10pNPfines-dr_DurnerFF_try6-004',
         'PresPlate-10%PF','PresPlate-15%PF','PresPlate-20%PF','PresPlate-30%PF']
TESTSETS = [['10% fines'   , 'galaa-10pfines-dr_DurnerFF_try6-006'],
            ['15% fines'   , 'galaa-15pfines-dr_DurnerFF_try7-002'],
            ['20% fines'   , 'galaa-20pfines-dr_DurnerFF_try6-002'],
            ['10% NP fines', 'galaa-10pNPfines-dr_DurnerFF_try6-004'],
            ['20% NP fines', 'galaa-20pNPfines-dr_DurnerFF_try5-002'],
            ['30% NP fines', 'galaa-30pNPfines-dr_DurnerFF_try6-002'],
           ]
TESTSETS = [['galaa-10pfines-dr_DurnerFF_try6-006',
             'galaa-15pfines-dr_DurnerFF_try7-002',
             'galaa-20pfines-dr_DurnerFF_try6-002',
             'galaa-10pNPfines-dr_DurnerFF_try6-004',
             'galaa-20pNPfines-dr_DurnerFF_try5-002',
             'galaa-30pNPfines-dr_DurnerFF_try6-002'],
            ['PresPlate-10%PF', 'galaa-10pfines-dr_DurnerFF_try6-006'],
            ['PresPlate-15%PF', 'galaa-15pfines-dr_DurnerFF_try7-002'],
            ['PresPlate-20%PF', 'galaa-20pfines-dr_DurnerFF_try6-002'],
           ]
#TESTSETS = [
#            ['PresPlate-10%PF', 'galaa-10pfines-dr_DurnerFF_try6-006'],
#           ]
SINGLE_PLOT = False

expdone = {}
expdone['10% fines'] = {
    'type' : 'Durner',
    'label': '10% PF - est',
    'ls': '-', #'r-', 'c-', 'g-'
    'color': 'b',
    'Ks': 5.82E-7,
    'n1': 29.497,
    'gamma1': -0.631,
    'n2': 1.187,
    'gamma2': -0.024,
    'w1': 0.165}
expdone['15% fines'] = {
    'type' : 'Durner',
    'label': '15% PF - est',
    'ls': '-', #'r-', 'c-', 'g-'
    'color': 'r',
    'Ks': 2.21E-7,
    'n1': 1.111,
    'gamma1': -27.403,
    'n2': 1.627,
    'gamma2': -0.003,
    'w1': 0.156}
expdone['20% fines'] = {
    'type' : 'vG',
    'label': '20% PF - est',
    'ls': '-', #'r-', 'c-', 'g-'
    'color': 'g',
    'Ks': 1.23E-7,
    'n': 1.171,
    'gamma': -0.039,}
expdone['10% NP fines'] = {
    'type' : 'Durner',
    'label': '10% NPF - est',
    'ls': '-', #'r-', 'c-', 'g-'
    'color': 'c',
    'Ks': 3.00E-8,
    'n1': 1.173,
    'gamma1': -0.535,
    'n2': 1.839,
    'gamma2': -0.022,
    'w1': 0.208}
expdone['20% NP fines'] = {
    'type' : 'Durner',
    'label': '20% NPF - est',
    'ls': '-', #'r-', 'c-', 'g-'
    'color': 'y',
    'Ks': 7.93E-9,
    'n1': 47.077,
    'gamma1': -0.646,
    'n2': 2.011,
    'gamma2': -0.012,
    'w1': 0.171}
expdone['30% NP fines'] = {
    'type' : 'Durner',
    'label': '30% NPF - est',
    'ls': '-', #'r-', 'c-', 'g-'
    'color': 'm',
    'Ks': 7.17E-9,
    'n1': 5.926,
    'gamma1': -0.0015,
    'n2': 1.321 ,
    'gamma2': -0.035,
    'w1': 0.229,}

expdone['20pNPF_DurnerFF_try4'] = {
    'type' : 'DurnerFF',
    'label': '20% NPF',
    'ls': '-', #'r-', 'c-', 'g-'
    'color': 'k',
    'gamma2'   : -0.01198936,
    'gamma1'   : -0.6459756,
    'Ks'       : 7.91674655e-07,
    'w1'       : 0.17095876,
    'n1'       : 47.21602695,
    'n2'       : 2.01094139,
    #hiorig =  [-800.,  -184.69972153,  -69.04066976,  -84.88393601,  -25.55339634, -1.]
    'hiorig' : [-800., -184.69972153, -84.88393601, -25.55339634, -1.],
    #kiorig = [  2.55592585e-177, 8.37197116e-084, 2.45017122e-077,   2.41964065e-004,
    #            1.89611978e-003,   4.73289371e-001]
    'kiorig' : [  4.28434076e-177,   1.54026912e-083,   7.29970690e-005,
                2.33028423e-003,   7.03968739e-001],

    #hinew = [-1198.8, -885.5, -498.7, -124.9, -55.6, -14.0, -1.6, -0.01]
    'hinew' : [-1198.8, -885.5, -498.7, -124.9, -110, -90, -70, -55.6, -14.0, -1.6, -0.01],
}

expdone['galaa-10pfines-dr_DurnerFF_try6-006'] = {
    'type' : 'DurnerFF',
    'label': '10% PF',
    'ls': '-',
    'color': 'k',
    'gamma2'   : -0.024,
    'gamma1'   : -0.631,
    'Ks'       :  5.82000000e-05,
    'w1'       : 0.165,
    'n1'       : 29.497,
    'n2'       : 1.187,
    #hiorig =  [-800.,  -184.69972153,  -69.04066976,  -84.88393601,  -25.55339634, -1.]
    'hiorig' : [-800.,         -353.78609229, -346.89955171,
                -194.86321058, -150.43565068, -144.52534243,
                -29.77738806 ,  -1.],
    #kiorig = [  2.55592585e-177, 8.37197116e-084, 2.45017122e-077,   2.41964065e-004,
    #            1.89611978e-003,   4.73289371e-001]
    'kiorig' : [0.00031817,  0.00101114,  0.00101148,
                0.001011484,  0.00637593,  0.01008568,
                0.01740192,  0.28190157],
}
expdone['galaa-15pfines-dr_DurnerFF_try7-002'] = {
    'type' : 'DurnerFF',
    'label': '15% PF',
    'ls': ':', #'r-', 'c-', 'g-'
    'color': 'k',
    'gamma2'   : -0.00444914,
    'gamma1'   : -7.30294145,
    'Ks'       :  5.06861494e-05,
    'w1'       : 0.20518757,
    'n1'       : 1.38017945,
    'n2'       : 1.49017759,
    #hiorig =  [-800.,  -184.69972153,  -69.04066976,  -84.88393601,  -25.55339634, -1.]
    'hiorig' : [-800.,     -161.95143678, -19.48271583,  -1.],
    #kiorig = [  2.55592585e-177, 8.37197116e-084, 2.45017122e-077,   2.41964065e-004,
    #            1.89611978e-003,   4.73289371e-001]
    'kiorig' : [1.05750234e-04,   5.81426318e-03,
                7.81000648e-02,   2.84713766e-01],

}
expdone['galaa-20pfines-dr_DurnerFF_try6-002'] = {
    'type' : 'DurnerFF',
    'label': '20% PF',
    'ls': '--', #'r-', 'c-', 'g-'
    'color': 'k',
    'gamma2'   : -0.04214436,
    'gamma1'   : -0.05749302,
    'Ks'       :  5.26936643e-05,
    'w1'       : 0.01242649,
    'n1'       : 1.19134601,
    'n2'       : 1.17209924,
    #hiorig =  [-800.,  -184.69972153,  -69.04066976,  -84.88393601,  -25.55339634, -1.]
    'hiorig' : [-800.,         -234.09078708,
                -37.81858234,   -1.],
    #kiorig = [  2.55592585e-177, 8.37197116e-084, 2.45017122e-077,   2.41964065e-004,
    #            1.89611978e-003,   4.73289371e-001]
    'kiorig' : [0.00326819,  0.00373306,  0.14241173,
                0.65473989],
}
expdone['galaa-10pNPfines-dr_DurnerFF_try6-004'] = {
    'type' : 'DurnerFF',
    'label': '10% NPF',
    'ls': '-', #'r-', 'c-', 'g-'
    'color': 'silver',
    'gamma2'   : -0.03936794,
    'gamma1'   : -0.58575879,
    'Ks'       :  8.80728887e-06,
    'w1'       : 0.22588523,
    'n1'       : 1.29315969,
    'n2'       : 1.8561588,
    #hiorig =  [-800.,  -184.69972153,  -69.04066976,  -84.88393601,  -25.55339634, -1.]
    'hiorig' : [-800.,    -347.49526585, -151.299384,
                -150.74744048,  -35.49824755,   -1.],
    #kiorig = [  2.55592585e-177, 8.37197116e-084, 2.45017122e-077,   2.41964065e-004,
    #            1.89611978e-003,   4.73289371e-001]
    'kiorig' : [0.02422905,  0.05657715,  0.05953699,
                0.05965733,  0.42722635,  0.87357708],
}
expdone['galaa-20pNPfines-dr_DurnerFF_try5-002'] = {
    'type' : 'DurnerFF',
    'label': '20% NPF',
    'ls': '--', #'r-', 'c-', 'g-'
    'color': 'silver',
    'gamma2'   : -0.01493587,
    'gamma1'   : -0.63955538,
    'Ks'       :  3.75972604e-05,
    'w1'       : 0.17469983,
    'n1'       : 46.02351962,
    'n2'       : 2.01108216,
    #hiorig =  [-800.,  -184.69972153,  -69.04066976,  -84.88393601,  -25.55339634, -1.]
    'hiorig' : [-800.,         -184.69972153,
                -30.01083991,   -1.],
    #kiorig = [  2.55592585e-177, 8.37197116e-084, 2.45017122e-077,   2.41964065e-004,
    #            1.89611978e-003,   4.73289371e-001]
    'kiorig' : [2.82850568e-04,   3.80381402e-03,
                4.03184581e-02,   6.49784135e-01],
}
expdone['galaa-30pNPfines-dr_DurnerFF_try6-002'] = {
    'type' : 'DurnerFF',
    'label': '30% NPF',
    'ls': '-.', #'r-', 'c-', 'g-'
    'color': 'silver',
    'gamma2'   : -0.06495382,
    'gamma1'   : -0.03686251,
    'Ks'       : 1.79724123e-06,
    'w1'       : 0.03948866,
    'n1'       : 2.30222211,
    'n2'       : 1.56687515,
    #hiorig =  [-800.,  -184.69972153,  -69.04066976,  -84.88393601,  -25.55339634, -1.]
    'hiorig' : [-800.,         -273.57063589,
                -60.92499156,   -1.],
    #kiorig = [  2.55592585e-177, 8.37197116e-084, 2.45017122e-077,   2.41964065e-004,
    #            1.89611978e-003,   4.73289371e-001]
    'kiorig' : [0.00102784,  0.04231832,
                0.05591445,  0.75044259],
}
expdone['PresPlate-20%PF'] = {
    'type': 'SWRCpoints',
    'label': 'PP 20% PF',
    'ls': 'x', #'r-', 'c-', 'g-'
    'color': 'b',
    'hiorig': [10,	31.6227766017,	50.1187233627,	63.095734448,
               100,	316.2277660168,	1000,	5011.8723362727,
               15848.9319246111],
    'uiorig': [1.0835354428,	1.0448574607,	1.0147731519,	0.9769370225,
               0.9306961765,	1.0433260192,	0.9170641716,	0.7474757274,
               0.0785668882],
}
expdone['PresPlate-10%PF'] = {
    'type': 'SWRCpoints',
    'label': 'PP 10% PF',
    'ls': 'x', #'r-', 'c-', 'g-'
    'color': 'b',
    'hiorig': [10,	31.6227766017,	50.1187233627,	63.095734448,
               100,	316.2277660168,	1000,	5011.8723362727,
               15848.9319246111],
    'uiorig': [1.0088572071,	0.9146880268,	0.8354976821,	0.7939658827,
               0.7675707206,	0.7744358486,	0.6917030211,	0.4387925351,
               0.0419470521],
}
expdone['PresPlate-30%PF'] = {
    'type': 'SWRCpoints',
    'label': 'PP 30% PF',
    'ls': 'x', #'r-', 'c-', 'g-'
    'color': 'b',
    'hiorig': [10,	31.6227766017,	50.1187233627,	63.095734448,
               100,	316.2277660168,	1000,	5011.8723362727,
               15848.9319246111],
    'uiorig': [1.2983478725,	1.2709426605,	1.2588595746,	1.2512955198,
               1.2391137192,	1.3273000344,	1.2220997871,	1.0082958729,
               0.1142061397],
}
expdone['PresPlate-15%PF'] = {
    'type': 'SWRCpoints',
    'label': 'PP 15% PF',
    'ls': 'x', #'r-', 'c-', 'g-'
    'color': 'b',
    'hiorig': [10,	31.6227766017,	50.1187233627,	63.095734448,
               100,	316.2277660168,	1000,	5011.8723362727,
               15848.9319246111],
    'uiorig': [1.1156882908,	1.0765678959,	1.058338383,	1.0443927097,
               1.024384112,	1.1062883131,	0.9636585738,	0.6834009334,
               0.0636839603],
}

print ('Available experiments:', expdone.keys())

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


SWRCS = {}
for TEST in TESTS:
    if not TEST in expdone:
        raise ValueError, 'Unknown test %s' % TEST
    if expdone[TEST]['type'] == 'vG':
        gamma   = expdone[TEST]['gamma']
        Ks      = expdone[TEST]['Ks']
        n       = expdone[TEST]['n']

        SWRCS[TEST] = SC_vanGenuchten(n=n, gamma=gamma)
    elif expdone[TEST]['type'] == 'SWRCpoints':
        SWRCS[TEST] = {'hiorig': np.asarray(expdone[TEST]['hiorig']),
                       'uiorig': np.asarray(expdone[TEST]['uiorig']),
                      }
    elif expdone[TEST]['type'] == 'Durner':
        gamma2   = expdone[TEST]['gamma2']
        gamma1   = expdone[TEST]['gamma1']
        Ks       = expdone[TEST]['Ks']
        w1       = expdone[TEST]['w1']
        n1       = expdone[TEST]['n1']
        n2       = expdone[TEST]['n2']
        if 'a1' in expdone[TEST]:
            a1       = expdone[TEST]['a1']
            a2       = expdone[TEST]['a2']
        else:
            a1 = 0.5
            a2 = 0.5

        SWRCS[TEST] = SC_Durner(n1=n1, gamma1=gamma1, n2=n2, gamma2=gamma2,
                                  w1=w1, a1=a1, a2=a2)


    elif expdone[TEST]['type'] == 'DurnerFF':
        gamma2   = expdone[TEST]['gamma2']
        gamma1   = expdone[TEST]['gamma1']
        Ks       = expdone[TEST]['Ks']
        w1       = expdone[TEST]['w1']
        n1       = expdone[TEST]['n1']
        n2       = expdone[TEST]['n2']
        hiorig = expdone[TEST]['hiorig']
        kiorig = expdone[TEST]['kiorig']

        hiorig = np.asarray(hiorig)
        kiorig = np.asarray(kiorig)

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

        SWRCS[TEST] = SC_Durner_freeform(n1=n1, gamma1=gamma1, n2=n2, gamma2=gamma2,
                                  w1=w1, hi=hiorig, ki=kiorig)

        if REFINE:
            logh2logk = MonoCubicInterp(hnodesorig, kvalsorig)
            hinew  = expdone[TEST]['hinew']
            hinew = np.asarray(hinew)

            hnodesnew = np.empty(len(hinew)+2, float)
            kvalsnew = np.empty(len(hinew)+2, float)
            hnodesnew[1:-1] = -1*np.log(-hinew[:])
            hnodesnew[0] = -1*np.log(-hiorig[0]*1e3)
            hnodesnew[-1] = -1*np.log(-hmax)
            kvalsnew = logh2logk(hnodesnew)


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

            SWRCS[TEST+'_new'] = SC_Durner_freeform(n1=n1, gamma1=gamma1, n2=n2, gamma2=gamma2,
                                  w1=w1, hi=hinew, ki=kinew)
    else:
        raise ValueError, 'Unknown type %s' % expdone[TEST]['type']


fig = 1

hmax = 5.5e3
haxis = np.linspace(-1500., -0.0001, hmax)

if SINGLE_PLOT:
    for TEST in TESTS:
        SWRCorig = SWRCS[TEST]
        uaxorig = SWRCorig.h2u(haxis)
        kaxorig = SWRCorig.h2Kh(haxis, expdone[TEST]['Ks'])
        if REFINE:
            SWRC = SWRCS[TEST+'_new']
            uax = SWRC.h2u(haxis)
            kax = SWRC.h2Kh(haxis, expdone[TEST]['Ks'])

        pylab.figure(fig)
        fig +=1
        #pylab.ylim(ymax=-hinew[0])
        pylab.title("$S_e$ " + TEST)
        pylab.ylabel('Negative Head [cm]')
        pylab.yscale('log')
        pylab.xlabel('Saturation')
        pylab.plot(uaxorig, -haxis, 'k-', label="Orig")
        if REFINE:
            pylab.plot(uax, -haxis, 'b-', label="New")
            pylab.legend()


        pylab.figure(fig)
        fig +=1
        pylab.title("Hydraulic Conductivity $K(S_e)$ " + TEST)
        pylab.ylabel('Hydraulic Conductivity [cm/s]')
        pylab.yscale('log')
        pylab.xlabel('Saturation')
        pylab.ylim(ymin=1e-20)
        pylab.plot(uaxorig, kaxorig, 'k-', label="Orig")
        if REFINE:
            pylab.plot(uax, kax, 'b-', label="New")
            pylab.legend()

        pylab.figure(fig)
        fig +=1
        pylab.title("$K(h)$ " + TEST)
        pylab.ylabel('Hydraulic Conductivity [cm/s]')
        pylab.xlabel('Negative Head [cm]')
        pylab.yscale('log')
        pylab.xscale('log')
        pylab.ylim(ymin=1e-40)
        pylab.plot(-haxis, kaxorig, 'k-', label="Orig")
        if REFINE:
            pylab.plot(-haxis, kax, 'b-', label="New")
            pylab.plot(-hinew, kinew*expdone[TEST]['Ks'], 'kx', label="Grid")
            pylab.legend()

else:
    #combined plot to see changes
    for SET in TESTSETS:
        SETSWRC = {}
        for TEST in SET:
            if isinstance(SWRCS[TEST], dict):
                continue
            SETSWRC[TEST] = {}
            SETSWRC[TEST]['SWRC'] = SWRCS[TEST]
            SETSWRC[TEST]['uax']  = SWRCS[TEST].h2u(haxis)
            SETSWRC[TEST]['kax']  = SWRCS[TEST].h2Kh(haxis, expdone[TEST]['Ks'])


        pylab.figure(fig)
        fig +=1
        #pylab.ylim(ymax=-hinew[0])
        pylab.title("$S_e$ ")
        pylab.ylabel('Negative Head [cm]')
        pylab.yscale('log')
        pylab.ylim(ymin=1e-1, ymax=hmax)
        pylab.xlabel('Saturation')
        indst = 0
        for TEST in SET:
            if isinstance(SWRCS[TEST], dict):
                pylab.plot(SWRCS[TEST]['uiorig'],  SWRCS[TEST]['hiorig'],
                           expdone[TEST]['ls'],
                           color=expdone[TEST]['color'],
                           label=expdone[TEST]['label'])
            else:
                pylab.plot(SETSWRC[TEST]['uax'],  -haxis, expdone[TEST]['ls'],
                           color=expdone[TEST]['color'],
                           label=expdone[TEST]['label'])
            indst += 1
        pylab.legend(loc=3)


        pylab.figure(fig)
        fig +=1
        pylab.title("Hydraulic Conductivity $K(S_e)$ ")
        pylab.ylabel('Hydraulic Conductivity [cm/s]')
        pylab.yscale('log')
        pylab.xlabel('Saturation')
        pylab.ylim(ymin=1e-10)
        indst = 0
        for TEST in SET:
            if isinstance(SWRCS[TEST], dict):
                continue
            pylab.plot(SETSWRC[TEST]['uax'], SETSWRC[TEST]['kax'], expdone[TEST]['ls'],
                       color=expdone[TEST]['color'],
                       label=expdone[TEST]['label'])
            indst += 1
        pylab.legend(loc=2)

        pylab.figure(fig)
        fig +=1
        pylab.title("Hydraulic Conductivity $K(h)$ ")
        pylab.ylabel('Hydraulic Conductivity [cm/s]')
        pylab.xlabel('Negative Head [cm]')
        pylab.yscale('log')
        pylab.xscale('log')
        pylab.ylim(ymin=1e-10)
        pylab.xlim(xmin=1, xmax=hmax)
        indst = 0
        for TEST in SET:
            if isinstance(SWRCS[TEST], dict):
                continue
            pylab.plot(-haxis, SETSWRC[TEST]['kax'], expdone[TEST]['ls'],
                       color=expdone[TEST]['color'],
                       label=expdone[TEST]['label'])
            indst += 1
        pylab.legend()

pylab.show()

