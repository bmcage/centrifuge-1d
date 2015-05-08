from __future__ import division, print_function

import numpy as np
from scipy.integrate import simps

import sys

from scikits.odes.sundials.ida import IDA_RhsFunction
from ..shared.functions import y2x
from ..shared.solver import simulate_direct
from ..shared.show import show_results
from ..shared.consolidation import CON_GOMPERTZ

use_cons_water = True
CONVERT = False

def zida2e(z, model):
    #convert internal z to true void ratio
    # void ratio between  CON._e0 - CON._c and CON._e0
    lowerb = model.CON._e0 - model.CON._c
    upperb = model.CON._e0
    (first_idx, last_idx) = (model.first_idx, model.last_idx)

    return lowerb + (np.sin(z[first_idx:last_idx+1])+1)*(upperb-lowerb)/2

def zdotida2edot(z, zdot, model):
    #convert internal z to true void ratio
    # void ratio between  CON._e0 - CON._c and CON._e0
    lowerb = model.CON._e0 - model.CON._c
    upperb = model.CON._e0
    (first_idx, last_idx) = (model.first_idx, model.last_idx)

    return np.cos(z[first_idx:last_idx+1])*zdot[first_idx:last_idx+1]*(upperb-lowerb)/2

def e2zida(e, model, z):
    lowerb = model.CON._e0 - model.CON._c
    upperb = model.CON._e0
    (first_idx, last_idx) = (model.first_idx, model.last_idx)
    z[first_idx:last_idx+1] = np.arcsin(2*(e-lowerb)/(upperb-lowerb)-1)

class centrifuge_residual(IDA_RhsFunction):
    def __init__(self):
        IDA_RhsFunction.__init__(self)
        self.ecopy = None

    def evaluate(self, t, z, zdot, result, model):
        # F(t,h,dh/dt) = porosity * du/dh * dh/dt
        #                + d[K(h)*(dh/dr - omega^2/g * r)]/dr
        try:

            (rE, Lorig, porosityorig, eorig, fl2, fl1, gamma_w, gamma_s) = \
              (model.re, model.l0, model.porosity, model.e0, model.fl2, model.fl1,
                model.gamma_w, model.gamma_s)
            eorignumer = model.numerfact_e0 * eorig

            if fl1 > 0:
                print ('ERROR: porous filter at top not supported yet!')
                return 1

            (ell, L) = (z[model.wl_idx], z[model.L_idx])
            r0 = rE - fl2 - L

            rb_type = model.rb_type

            if (rb_type == 2):
                print('ERROR: rb_type is 2, but prescribed pressure at end sample not supported yet!')
                return 1
            if (rb_type == 3):
                print('ERROR: rb_type is 3, but overflow lip at end sample not supported yet!')
                return 1

            (first_idx, last_idx) = (model.first_idx, model.last_idx)
            #consolidation data
            CON = model.CON

            if self.ecopy is None:
                self.ecopy = np.empty(last_idx+1-first_idx, float)
            if CONVERT:
                self.ecopy[:] = zida2e(z, model)
            else:
                self.ecopy[:] = z[first_idx:last_idx+1]

            if CON.typeCON() == CON_GOMPERTZ:
                if np.any(self.ecopy <= CON._e0 - CON._c): # negative void ratio??
                    print ('ERROR: Too low void ratio found', self.ecopy, '<', CON._e0 - CON._c)
                    # added a regularization here. IMPORTANT
                    t1 = 1.001 - 0.0001*(CON._e0 - CON._c-self.ecopy[self.ecopy <= CON._e0 - CON._c])
                    t2 = 1.0001 - 0.0000001*(CON._e0 - CON._c-self.ecopy[self.ecopy <= CON._e0 - CON._c] )
                    self.ecopy[self.ecopy <= CON._e0 - CON._c] = \
                        np.max([np.max(t1),np.max(t2),1]) *  (CON._e0 - CON._c)
                    print ('Correcting, new value: ', self.ecopy)

            else:
                if np.any(self.ecopy <= 0): # negative void ratio??
                    print ('ERROR: Negative void ratio found', self.ecopy)
                    print ('input z   ', z)
                    print ('input zdot', zdot)
                    self.ecopy[self.ecopy <= 0] = 0

            if CONVERT:
                edot = zdotida2edot(z, zdot, model)
            else:
                edot =  zdot[first_idx:last_idx+1]

            e12  = (self.ecopy[1:] + self.ecopy[:-1]) / 2

            #constitutive laws
            Ks_e = CON.e2Ks(self.ecopy)
            Ks_e12 = CON.e2Ks(e12)

            if (t<1e-10):
                effstress_e = CON.e2sigmaprime(self.ecopy/model.numerfact_e0)
            else:
                effstress_e = CON.e2sigmaprime(self.ecopy)

            #the grid
            (y, y12, dy, alpha)  = (model.y, model.y12, model.dy, model.alpha)
            (ldc1, ldc2, ldc3) = (model.ldc1, model.ldc2, model.ldc3)

            dedy12 = (self.ecopy[1:] - self.ecopy[:-1]) / dy
            dedy = np.empty([model.inner_points+2,], dtype=float)
            dedy[0]    = ldc1[0] * self.ecopy[0] + ldc2[0] * self.ecopy[1] + ldc3[0]*self.ecopy[2]
            dedy[1:-1] = ldc1[1:-1]*self.ecopy[:-2] + ldc2[1:-1]*self.ecopy[1:-1] + ldc3[1:-1]*self.ecopy[2:]
            dedy[-1]   = ldc1[-1] * self.ecopy[-3] + ldc2[-1] * self.ecopy[-2] + ldc3[-1] * self.ecopy[-1]

            omega2g = model.find_omega2g(t)
            #small regularization to start with some effect
            minom = 0#1e-3
            minomega2g = minom*minom/model.g
            if omega2g < minomega2g:
                omega2g = minomega2g

            (delldt, dLdt) = (zdot[model.wl_idx], zdot[model.L_idx])

            #derivatives constitutive laws
            dKo1pe_de = CON.dKo1pede(self.ecopy)  #d (K/(1+e)) / de at e_i
            dsigp_de12  = CON.dsigpde(e12) #d sigma' / de at e_{i+1/2}

            #no filter at top, BC at top is a water level ell(t)
            load = model.excess_load_f(t)

            if load == 0:
                result[first_idx] = eorignumer - self.ecopy[0]  # at top, no consolidation!
            else:
                # at top, the effective stress is the load
                result[first_idx] = CON.sigmaprime2e(load) - self.ecopy[0]

            segment12 = (alpha[1:-1] + alpha[2:])/2
            Dpart = Ks_e12 * dsigp_de12 * dedy12

            result[first_idx+1:last_idx] = \
              segment12 * edot[1:-1] - segment12 *(1-y[1:-1])/L*dLdt*dedy[1:-1] \
             -segment12 * (gamma_s/gamma_w-1)*omega2g*(1+self.ecopy[1:-1]) \
                        * (rE-fl2*(1-y[1:-1])*L)/L * dKo1pe_de[1:-1] * dedy[1:-1] \
             -segment12 * (gamma_s/gamma_w-1)*omega2g*Ks_e[1:-1] \
             + (1+self.ecopy[1:-1])/gamma_w/L/L * ( Dpart[1:] - Dpart[:-1])

            #we want e to remain below max, so add an error if above
            if model.e0_overshoot_factor:
                bade = eorig < self.ecopy[1:-1]
                result[first_idx+1:last_idx][bade] += model.e0_overshoot_factor*(self.ecopy[1:-1][bade] - eorig) \
                                * np.sign(result[first_idx+1:last_idx][bade])

            #pressure at top
            pL = gamma_w*omega2g/2*( (rE-fl2-L-fl1)**2 - (rE-fl2-L-fl1*ell)**2 )
            #total stress over the sample
            totsig = np.empty(len(self.ecopy), float)
            totsig[0] = pL + load
            for i in range(1,len(self.ecopy)):
                totsig[i] = totsig[i-1] \
                    + omega2g*(gamma_s+gamma_w*e12[i-1])/(1+e12[i-1]) \
                             *(rE-fl2-(1-y12[i-1])*L) * L * dy[i-1]

            if rb_type == 3:
                return 1  # not yet programmed, see run of direct_sat_drain
            elif rb_type == 0:  #no outflow!
                dsigde = CON.dsigpde([self.ecopy[-1]],zeroval=-1e-20)
                if CON.typeCON() == CON_GOMPERTZ:
                    eatend = self.ecopy[-1]
                    if eatend>= CON._e0:
                        eatend=CON._e0

                    # problem, at edge we have e>=e0, and sig == 0, dsigde ==0!
                    # we estimate dsigdy instead
                    if (self.ecopy[-3] <= CON._e0 and self.ecopy[-2] <= CON._e0
                            and self.ecopy[-1] <= CON._e0):
                        #valid void ratio's
                        dsigdy   = ldc1[-1] * CON.e2sigmaprime([self.ecopy[-3]]) + ldc2[-1] * CON.e2sigmaprime([self.ecopy[-2]]) + ldc3[-1] * CON.e2sigmaprime([self.ecopy[-1]])
                    elif (self.ecopy[-1] >= CON._e0):
                        #invalid void ratio, we cannot know a corresponding dsigdy
                        #we return dsigdy as -dedy to force solver to correct and move
                        #to values below e0
                        dsigdy = - ((self.ecopy[-1]-self.ecopy[-2])/dy[-1])
                    else:
                        #numerical discretization based on the last two values, last will be < e0
                        dsigdy = (CON.e2sigmaprime([self.ecopy[-1]])-CON.e2sigmaprime([self.ecopy[-2]]))/dy[-1]

                    #boundary condition sets the value of dsigdy
                    result[last_idx]  = ( dsigdy
                            - (L
                               * (gamma_s-gamma_w)
                               / (1+eatend) *omega2g * (rE-fl2))
                                      )
                    if (self.ecopy[-1] >= CON._e0):
                        #we had step in impossible direction. penelize it!
                        result[last_idx] *= (1+self.ecopy[-1] - CON._e0)
#==============================================================================
#                     if dsigde == 0:
#                         # problem, at edge we have e>=e0, and sig == 0, dsigde ==0!
#                         # we estimate dsigdy instead
#                         if (self.ecopy[-3] <= CON._e0 and self.ecopy[-2] <= CON._e0 and self.ecopy[-1] <= CON._e0):
#                             #valid void ratio's
#                             dsigdy   = ldc1[-1] * CON.e2sigmaprime([self.ecopy[-3]]) + ldc2[-1] * CON.e2sigmaprime([self.ecopy[-2]]) + ldc3[-1] * CON.e2sigmaprime([self.ecopy[-1]])
#                         elif (self.ecopy[-1] >= CON._e0):
#                             #invalid void ratio, we cannot know a corresponding dsigdy
#                             #we return dsigdy as -dedy to force solver to correct
#                             dsigdy = - ((self.ecopy[-1]-self.ecopy[-2])/dy[-1])
#                         else:
#                             dsigdy = (CON.e2sigmaprime([self.ecopy[-1]])-CON.e2sigmaprime([self.ecopy[-2]]))/dy[-1]
#
#                         result[last_idx]  = ( dsigdy
#                                 - (L
#                                    * (gamma_s-gamma_w)
#                                    / (1+eatend) *omega2g * (rE-fl2))
#                                           )
#                         if (self.ecopy[-1] >= CON._e0):
#                             #we had step in impossilbe direction. penelize it!
#                             result[last_idx] *= (1+self.ecopy[-1] - CON._e0)
#                     else:
#                         result[last_idx]  = (1e-6+omega2g)*( dsigde
#                             * dedy[-1]
#                                 - (L
#                                    * (gamma_s-gamma_w)
#                                    / (1+self.ecopy[-1]) *omega2g * (rE-fl2))
#                                           )
#                         dsigdy   = ldc1[-1] * CON.e2sigmaprime([self.ecopy[-3]]) + ldc2[-1] * CON.e2sigmaprime([self.ecopy[-2]]) + ldc3[-1] * CON.e2sigmaprime([self.ecopy[-1]])
#                         #dsigdy = (CON.e2sigmaprime([self.ecopy[-1]])-CON.e2sigmaprime([self.ecopy[-2]]))/dy[-1]
#                         if (self.ecopy[-1] >= CON._e0):
#                             #invalid void ratio, we cannot know a corresponding dsigdy
#                             #we return dsigdy as -dedy to force solver to correct
#                             dsigdy = - ((self.ecopy[-1]-self.ecopy[-2])/dy[-1])
#                         result[last_idx]  = ( dsigdy
#                                 - (L
#                                    * (gamma_s-gamma_w)
#                                    / (1+self.ecopy[-1]) *omega2g * (rE-fl2))
#                                          )
#==============================================================================
                else:
                    result[last_idx]  = (1e-6+omega2g)*( dsigde
                        * dedy[-1]
                            - (L
                               * (gamma_s-gamma_w)
                               / (1+self.ecopy[-1]) *omega2g * (rE-fl2))
                                      )
                    #dsigdy   = ldc1[-1] * CON.e2sigmaprime([self.ecopy[-3]]) + ldc2[-1] * CON.e2sigmaprime([self.ecopy[-2]]) + ldc3[-1] * CON.e2sigmaprime([self.ecopy[-1]])
                    #dsigdy = (CON.e2sigmaprime([self.ecopy[-1]])-CON.e2sigmaprime([self.ecopy[-2]]))/dy[-1]
                    #result[last_idx]  = ( dsigdy
                    #        - (L
                    #           * (gamma_s-gamma_w)
                    #           / (1+self.ecopy[-1]) *omega2g * (rE-fl2))
                    #                  )

#==============================================================================
#                 if omega2g < 0.5e-7:
#                     # too small rotation, no change yet!
#                     result[last_idx] = dedy[-1]
#                 elif omega2g < 1e-7:
#                     #regularize start up
#                     result[last_idx] = dedy[-1] + (omega2g-0.5e-7)/0.5e-7 * (result[last_idx]- dedy[-1])
#==============================================================================
            elif rb_type == 1: #free outflow
                #effective stress should be equal to total stress at BC
                #print ('Test free outflow' , t, effstress_e[-1], totsig[-1], CON.sigmaprime2e(totsig[-1]), np.power(omega2g*model.g,.5))
                result[last_idx]  = self.ecopy[-1] - CON.sigmaprime2e(totsig[-1])
                #result[last_idx]  = effstress_e[-1] - totsig[-1]
            elif rb_type == 2:
                return 1  #not yet programmed
            else:
                raise NotImplementedError('rb_type has to be 0 - 3')

            dsminsprimedy_left  = ldc1[0] * (totsig[0]-effstress_e[0]) \
                                + ldc2[0] * (totsig[1]-effstress_e[1]) \
                                + ldc3[0] * (totsig[2]-effstress_e[2])
            dsminsprimedy_right = ldc1[-1] * (totsig[-3]-effstress_e[-3]) \
                                + ldc2[-1] * (totsig[-2]-effstress_e[-2]) \
                                + ldc3[-1] * (totsig[-1]-effstress_e[-1])

            #length sample
            int_1_over_1_plus_e = simps(1/( 1+self.ecopy[:] ), x=y)
            #computed L as an algebraic equation (see solve, algvars_idx = [model.L_idx] ):
            result[model.L_idx] = L - Lorig/(1+eorignumer)/int_1_over_1_plus_e

            #mass flowing out. This is expressed in length as ell and L, so
            # we assume an area equal to area of the sample tube. This can
            # be corrected in postprocessing for gf_mo, not here !!!
            if rb_type == 0:
                #Mout must be fixed, so derivative should be zero
                result[model.mass_out_idx] = zdot[model.mass_out_idx]
            else:
                qN = -Ks_e[-1]/gamma_w * (1/L*dsminsprimedy_right \
                                          - gamma_w*omega2g*(rE-fl2))
                result[model.mass_out_idx] = zdot[model.mass_out_idx] - qN


            #water level on top
            ## NOTE: we alternatively can use an algebraic equation of
            ##       conservation of water !!
            if use_cons_water:
                (wm, wm_in_tube) = model.measurements.store_calc_wm_e(
                        self.ecopy,
                        model.y, z[model.L_idx], z[model.wl_idx],
                        z[model.mass_out_idx],
                        model.fl2, model.fp2, store=False)
                result[model.wl_idx] = wm - model.wm0
            else:
                q0 = -Ks_e[0]/gamma_w * (1/L*dsminsprimedy_left \
                                      - gamma_w*omega2g*(rE-fl2-L))
                result[model.wl_idx] = zdot[model.wl_idx] + q0

            #mass flowing in: not used at the moment in this model
            result[model.mass_in_idx] = zdot[model.mass_in_idx]

#            print ('test', t, result,[z[last_idx-2],z[last_idx-1], z[last_idx]],
#                        zdot[last_idx],
#                        [result[last_idx-2],result[last_idx-1],result[last_idx]])
#            #raw_input('cont')

            verbosity = model.verbosity
            if verbosity > 2:
                if verbosity > 3:
                    print('t: %10.6g' % t, 'wl: %8.6g' % z[model.wl_idx],
                          'L: %8.6g' % z[model.L_idx])
        except:
            import traceback

            exc_type, exc_value, exc_traceback = sys.exc_info()
            print ("*** print_exception:")
            traceback.print_exception(exc_type, exc_value, exc_traceback,
                                      limit=2, file=sys.stdout)
            return -1
#==============================================================================
#         print ('input z   ', z)
#         print ('input zdot', zdot)
#         print ('residual  ', result)
#         raw_input('TEST')
#==============================================================================
        return 0

RESIDUAL_FN = centrifuge_residual()
def sat_cons_root(t, z, zdot, result, model):
    #print ('cur root val', z[model.wl_idx])
    result[0] = z[model.wl_idx]

def sat_cons_root_cont(model, t, z):
    return False

def on_measurement(t, z, model, measurements):

    L  = measurements.store_calc_measurement('L', z[model.L_idx])
    wl = measurements.store_calc_measurement('wl', z[model.wl_idx])
    measurements.store_calc_measurement('MI', z[model.wl_idx])
    MO = measurements.store_calc_measurement('MO', z[model.mass_out_idx])

    #compute the x values from 0 to L corresponding to grid y values:
    x  = y2x(model.y, 0, L)

    if CONVERT:
        e = zida2e(z, model)
    else:
        e = z[model.first_idx: model.last_idx+1]

    (Ks_e, effstress_e) = measurements.store_calc_e(x, e, model.CON)

    (WM, WM_in_tube) = \
      measurements.store_calc_wm_e(e, model.y, L, wl, MO,
                                 model.fl2, model.fp2)

    if measurements.calc_measurement_p('GC'):
        measurements.store_calc_gc_e(e, model.y, model.dy, L,
                                     model.density, model.density_s)

    omega2g = model.find_omega2g(t)

    if measurements.calc_measurement_p('gF_MT'):
        measurements.store_calc_gf_mt_e(omega2g, e, model.y, model.dy, L,
                                      wl, model.re,
                                      model.fl2, model.fp2,
                                      model.density, model.density_s,
                                      model.tube_crosssectional_area)

    if measurements.calc_measurement_p('gF_MO'):
        measurements.store_calc_gf_mo(omega2g, MO,
                                      model.mo_gc_calibration_curve,
                                      model.density,
                                      model.tube_crosssectional_area)

def initialize_z0(z0, model):
    z0[model.first_idx:model.last_idx+1] = model.numerfact_e0 * model.e0

    mass_in = mass_out = 0.0

    z0[model.L_idx] = model.l0
    z0[model.wl_idx] = model.wl0
    z0[model.mass_in_idx]  = mass_in
    z0[model.mass_out_idx] = mass_out

    # assign value to WM0 and wm0 (wm0 is needed for mass balance)

    (wm0, wm_in_tube0) = \
      model.measurements.store_calc_wm_e(z0[model.first_idx:model.last_idx+1],
            model.y, model.l0, model.wl0, mass_out,
            model.fl2, model.fp2, store=False)

    model.wm0 = wm0

    if CONVERT:
        e2zida(z0[model.first_idx:model.last_idx+1], model, z0)


def initialize_zp0(zp0, z0, model):
    (first_idx, last_idx) = (model.first_idx, model.last_idx)
    #e = z0[first_idx:last_idx+1]

    t = 0.0
    omega2g = model.find_omega2g(t)

    if omega2g == 0:
        #nothing should change yet
        zp0[first_idx:last_idx+1][:] = 0
        zp0[model.L_idx] = 0
        zp0[model.wl_idx] = 0
        zp0[model.mass_in_idx] = 0
        zp0[model.mass_out_idx] = 0
    else:
        print('Not supported to start with a rotation non zero')
        exit(1)

def solve(model, measurements):
    measurements.reset_calc_measurements()

    # Initialization. L is via an algebraic equation based on e values
    if model.rb_type == 0:
        algvars_idx = [model.last_idx, model.L_idx]
    elif model.rb_type == 1:
        algvars_idx = [model.last_idx, model.L_idx]
    else:
        print('Not supported yet')
        sys.exit(1)
    if use_cons_water:
        algvars_idx += [model.wl_idx]

    atol_backup        = model.atol # backup value
    if type(atol_backup) == list:
        atol = np.asarray(atol_backup, dtype=float)
    else:
        atol = atol_backup * np.ones([model.z_size,], dtype=float)
    atol[model.L_idx] = model.L_atol
    model.atol         = atol

    if model.estimate_zp0:
        zp0_init = initialize_zp0
    else:
        zp0_init = None

    # Computation
    (flag, t, z, i) = simulate_direct(initialize_z0, model, model.measurements,
                                      RESIDUAL_FN,
                                      root_fn = sat_cons_root, nr_rootfns=1,
                                      continue_on_root_found = sat_cons_root_cont,
                                      truncate_results_on_stop = True,
                                      initialize_zp0=zp0_init,
                                      algvars_idx=algvars_idx,
                                    #  exclude_algvar_from_error=True,
                                      on_measurement=on_measurement)

    # Restore modified values
    if model.rb_type == 3:
        model.atol = atol_backup

    return flag

solve_direct = solve  # to remove warning about 'solve_direct' not specified

import cProfile  #, pstats
PROFILE = False

def run(model):
    if PROFILE:
        pr = cProfile.Profile()
        pr.enable()

    show_results(model.experiment_info, model=model)

    if PROFILE:
        pr.disable()
        pr.print_stats('cumtime')

def dry_run(model):
    show_results(model.experiment_info, model=model)
