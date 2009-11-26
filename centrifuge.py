#!/usr/bin env python

# Copyright (C) 2008  B. Malengier
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

import sys, os
import numpy as N
import math
from scipy import * 
from scipy.integrate import odeint,simps
import Gnuplot
import pylab as pl
import copy


#-------------------------------------------------------------------------
#
# local modules
#
#-------------------------------------------------------------------------
import UtilsBm
import GridUtils
from diffusion import *
from LinDiff import *


#-------------------
#
# module constant
#
#-------------------

g = 9.80665 #standard gravity, m/s^2
rho_f = 998.2071 # density at 20 degree Celcius, kg / m^3

class DataProblem():
    '''general class to hold data of a problem
    '''
    def __init__(self, options):
        ''' init the problem with options holding everything
        '''
        self.options = options
        self.info = '\n'
    
    def data(self):
        ''' return a string with the data'''
        return self.info
    def end(self):
        return self.options.endTime
    def nrsteps(self):
        return self.options.nrsteps
    def outputtime(self):
        return self.options.outputtime


class DataCentrifuge(DataProblem):
    '''data for a centrifuge problem
    '''
    def __init__(self, options):
        DataProblem.__init__(self, options)
        #check if all is present that is needed:
        self.omega_val = eval(self.options.omega)
        self.g_val = eval(self.options.g)
        
        self.info += """
Centrifuge flow with data:
    Axis at 0, with rotational speed %(omega)s, container at x=%(leftpos)f.
    In the container, sample starts at x=%(leftposL)f, with length %(L)f, 
     left porous stone of thickness %(th1)f, right of thickness %(th2)f. 
     Porous stone has internal permeability %(kp)g and porosity %(np)f. 
    The sample has porosity %(n)f, internal permeability %(k)g and 
     compressibility %(alpha)g.
""" % {'omega': self.options.omega, 'leftpos': self.leftpos(), 'L': self.L(),
            'leftposL': self.leftposL(), 'th1': self.th_p1(), 
            'th2': self.th_p2(), 'kp': self.k_p(), 'np': self.n_p(),
            'n': self.n0(), 'k': self.k(), 'alpha': self.alpha(),
            }
        
        self.BC()
        #only infiltrating free water column infiltration, so dirichlet at
        # boundary
        assert (self.BC() == 'dirdir' or self.BC() == 'dirhomdir' or
                self.BC() == 'freedirdir' or self.BC() == 'freedirhomdir' 
                )
        if self.BC()[:4] == 'free':
            # left is a free infiltrating water column
            self.bctype = 'free'
            leftbc = 'free inflow of water column, s(0) = %f' % self.options.s0
        elif self.BC()[:4] == 'func':
            # left is a prescribed infiltrating water column over time
            self.bctype = 'func'
            leftbc = 'prescribed function of height water column'
        else:
            self.bctype = ''
            leftbc = 'fixed prescribed pressure p1 = %f' % self.options.p1
        #only infiltrating water column implemented
        assert (self.pcalc() == 'h')
        if self.pcalc() == 'h':
            self.p1_val = None
            self.s0_val = self.options.s0
            self.h2_val = self.options.h2
            if self.h2_val is not None :
                self.p2_val = rho_f * g * self.h2_val
            else:
                self.p2_val = self.options.p2
                self.h2_val = self.p2_val / rho_f / g
        else :
            self.p1_val = self.options.p1
            self.p2_val = self.options.p2
            self.s0_val = None
        print self.p1_val,  self.s0_val,  self.p2_val,  self.h2_val
        self.info += """    Boundary conditions are: Left %(leftbc)s,
      Right %(rightbc)s
""" % {'leftbc': leftbc, 'rightbc': ' TODO '}
        self.info += """    The sample is at temperature %(temp)s, with water viscosity 
      %(vis)g, and water density %(dens)g.
""" % {'temp': self.temp(), 'vis': self.mu(), 'dens': self.rho_f()}

    def inittype(self):
        return self.options.inittype
    def k(self):
        return self.options.k
    def omega(self, t):
        return self.omega_val(t)
    def alpha(self):
        return self.options.alpha
    def n0(self):
        return self.options.n0
    def L(self):
        return self.options.L
    def leftposL(self):
        return self.options.leftposL
    def leftpos(self):
        return self.options.leftpos
    def g(self,t):
        return self.g_val(t)
    def th_p1(self):
        return self.options.th_p1
    def th_p2(self):
        return self.options.th_p2
    def k_p(self):
        return self.options.k_p
    def n_p(self):
        return self.options.n_p
    def BC(self):
        return self.options.BC
    def BCtype(self):
        return self.bctype
    def pcalc(self):
        return self.options.pcalc
    def p1(self):
        return self.p1_val
    def p2(self):
        return self.p2_val
    def s0(self):
        return self.s0_val
    def h2(self):
        return self.h2_val
    def mu(self):
        return self.options.mu
    def rho_f(self):
        return self.options.rho_f
    def temp(self):
        return self.options.temp
    def grav(self):
        return self.options.grav
    def alphagen(self):
        return self.options.alphagen
    def ngen(self):
        return self.options.ngen
    def theta_r(self):
        return self.options.theta_r
    def theta_s(self):
        return self.options.theta_s
    def alphagen_p(self):
        return self.options.alphagen_p
    def ngen_p(self):
        return self.options.ngen_p
    def theta_r_p(self):
        return self.options.theta_r_p
    def theta_s_p(self):
        return self.options.theta_s_p
    def eps(self):
        return self.options.eps

#-------------------------------------------------------------------------
#
# PDE solvers for centrifuge
#
#-------------------------------------------------------------------------

class CentrifugeFlow(Lin1DDiffusion):
    '''implement 1 D linear diffusion u_t = (D(x) u_x)_x for centrifuge
        input:  grid: the grid between [0,1]
                uvalinit : initial values over the grid
                data : Data of the sample
        Solution method: ....
    '''
    def __init__(self, grid, uvalinit, data, dir='test'):
        diff = Diffusion_const(1.)
        self.data = data
        #changable parameters
        self.k = self.data.k()
        self.alpha = self.data.alpha()
        Lin1DDiffusion.__init__(self, grid, uvalinit, diff, dir, data.L(), 
                                data.leftposL())

        self.lsodaprecr = 1e-5
        self.lsodapreca = 1e-2
        print '\n precision relative = %g, absolute = %g\n' % (self.lsodaprecr,
                                            self.lsodapreca)

        #override things to work with moving boundary
        if self.data.BCtype() == 'free':
            uvalinit_fv_real = zeros(len(self.uvalinit_fv)+1)
            uvalinit_fv_real[0] = self.data.s0()
            uvalinit_fv_real[1:] = self.uvalinit_fv[:]
            self.uvalinit_fv = uvalinit_fv_real
        
        #set time data for output
        self.timeexp = []
        self.nrexp = 0
        timestep  = self.data.end() / self.data.nrsteps()
        self.timeval   = []
        self.sol = []
        for i in arange(self.data.nrsteps()+1):
            self.timeval.append(i*timestep)
        self.timeval = sorted(UtilsBm.merge_nodup(self.timeval, self.timeexp))
        self.timeval = array(self.timeval)  #convert to numpy array
        
        
        #set BC data for infiltrating water column. Possibilities:
        #  1/free infiltration, height water column follows from infiltration
        #     ==> dirdir boundary with changing Dirichlet condition at inflow
        #  2/prescribed water height (for inverse problems)
        #     ==> neudir boundary with prescribed pressure change at inflow
        #         following from the water column decrease.
        
        if data.BC() == 'dirdir' or data.BC() == 'dirhomdir' \
                or data.BC() == 'freedirdir' or data.BC() == 'freedirhomdir':
            self.bc = 'dir'
            if data.BC() == 'dirhomdir' or data.BC() == 'freedirhomdir':
                assert self.data.p2() == 0.

##        elif bctype == 'homdirdir':
##            self.bc = 'dir'
##            self.bcpar = [0., self.data.p2()]
##        elif bctype == 'neuneu':
##            self.bc = 'neu'
##            self.bcpar = [bcpar[0], bcpar[1]]
##        elif bctype == 'homneuneu':
##            self.bc = 'neu'
##            self.bcpar = [0., bcpar[1]]
##        elif bctype == 'neuhomneu':
##            self.bc = 'neu'
##            self.bcpar = [bcpar[0], 0.]
##        elif bctype == 'homrobhomrob':  # du/dx = par *u 
##            self.bc = 'rob'
##            self.bcpar = [bcpar[0], 0., bcpar[1], 0.]
        else: print 'ERROR: Unknown border condition', data.BC();sys.exit()
        #override how bcval works
        self.bcvalparam = True

        #gradient at left and right boundary based on time and concentration in
        # cell adjacent to boundary
        if self.bc == 'neu':
            self.gradL = lambda t,u : self.bcpar[0]
            self.gradR = lambda t,u : self.bcpar[1]
        elif self.bc == 'rob':
            self.gradL = lambda t,u : self.bcpar[0] * u + self.bcpar[1]
            self.gradR = lambda t,u : self.bcpar[2] * u + self.bcpar[3]
        elif self.bc == 'dir':
            if self.data.pcalc() == 'h': 
                #deduce pressure at inflow from water table, omega can be variable!
                self.presL = lambda t: self.calc_presL(s=self.data.s0(),t=t)
            else:
                #given pressure, no change over time
                self.presL = lambda t: self.data.p1()
            if self.data.BCtype() == 'free':
                self.gradL = lambda t,s,u: (u-self.calc_presL(s,t))/self.Deltamin[0]
            else:
                self.gradL = lambda t,u : (u - self.presL(t))/self.Deltamin[0]
            self.gradR = lambda t,u : (self.data.p2() - u)/self.Deltamin[-1]
        
        # set possible parameters present, dictionary with get method, set
        #  method and method to return number of par for this type.
        #  Type is the dictionary key
        self.param_type = {
                'k' : [self.get_hydrcond, self.set_hydrcond
                                , self.return_one],
                'alpha' : [self.get_compr, self.set_compr
                                , self.return_one],
            }
        # new object, so changed parameters
        self.changedparam = True

    def return_one(self):
        return 1
    #methods to set and retrieve par
    def get_hydrcond(self, number) :
        return self.k
    def set_hydrcond(self, number, val) :
        self.k = val
    def get_compr(self, number) :
        return self.alpha
    def set_compr(self, number, val) :
        self.alpha = val
    def get_param(self, type, number) :
        return self.param_type[type][0](number)
    def set_param(self, type, number, val) :
        '''
        Inverse methods may only use this method, as it controls the
        changedparam attribute
        '''
        if self.param_type[type][0](number)!= val :
            self.changedparam = True
            self.param_type[type][1](number, val)
    def nrparam(self, type) :
        return self.param_type[type][2]()

    def calc_presL(self, s = 0., t=0.):
        """
        Calculate the pressure to the left given s, the position of the water
         table relative to begin of container and t, time of calculation
        We have a contribution of pressure of pure water, and a contribution
         of the water in the porous stone
        """
        contr = 0.
        rho_om2 = self.data.rho_f()* (self.data.omega(t) ** 2) / 2.
        n_p = self.data.n_p()
        pos_s = self.data.leftpos() + s
        posL = self.data.leftposL()
        #position of left position of porous stone.
        pos_pstL = posL - self.data.th_p1()
        #we do not include under pressure!
        assert (pos_s < posL)
        if pos_s < pos_pstL:
            #we have a contribution of pure water
            contr += pos_pstL **2 - pos_s **2
            #and of the porous stone
            contr += n_p * (posL** 2 - pos_pstL**2)
        else: #water in porous stone
            contr += n_p * (posL** 2 - pos_s**2)
        return contr * rho_om2

    def _solve(self, init, timeval):
        if self.solmeth == 'conservative_fv' :
             return odeint(self.f_cons_fv, copy.deepcopy(init),timeval, \
                            rtol=self.lsodaprecr, atol=self.lsodapreca, printmessg = 1, \
                            ml=1, mu=1)
        else :
            print 'ERROR : unknown solution method to solve ODE system'
            sys.exit()

    def f_cons_fv(self,u,t) :
        ''' 
            Overwrite standard lin1ddiffusion method
            right handed side of 1st - Order ODE system, u[0]=s, u[1:]=p
                s_t = -k/mu * (omega^2 leftposL rho_f + \partial_x p(leftposL))
                u_t = (\partial_x^2 p - rho_f omega^2)k/mu/alpha/(1-n)*rho_f*g*exp(alpha p)
                conservative form: 
                    u_t = factor * (flux_right - flux_left) + source
            We work with cell data, all cells are full cells
        '''
        np   = len(self.grid) #number of edges
        #it happens that the odeint is stopped and restarted, we need to have correct time:
        realtime = self.shifttime + t
        #set all inner points, uvalinit also np-1 large + edge movement
        u_t    = zeros(np, float)
        # an adjoint problem internally works with an inverted time
        # calculate the real time.
        #time = t
        #if self.adjoint :
        #    time = self.end - t

        #fluxes on all edges (derived from cell data)
        fluxmin = zeros(np,float)  #dp/dx
        fluxmin[1:-1] = (u[2:]-u[1:-1])/     \
                        ((self.Deltamin[:-1]+self.Deltamin[1:])/2.)
        #border flux left:
        fluxleftBC = self.gradL(realtime,u[0], u[1])
        #print 'fluxleftBC', fluxleftBC, t
        #normal to the left is -1!
        fluxmin[0] = fluxleftBC
        fluxrightBC = self.gradR(realtime, u[-1])
        fluxmin[-1] = fluxrightBC

        t1     = self.length * self.length
        t2     = 1./t1
        source = self.data.rho_f() * (self.data.omega(realtime)**2)
        factor = self.data.k()/self.data.mu()*self.data.rho_f() \
                    * self.data.grav()/self.data.alpha()/(1.-self.data.n0())\
                    * exp(self.data.alpha() * u[1:])
        u_t[1:] = (t2 * (fluxmin[1:]- fluxmin[:-1]) \
                        / (self.Deltamin[:]) - source) * factor[:]
        #speed of water table, diff(s(t),t)
        u_t[0] = -self.data.k()/self.data.mu() *  \
                    (( -1.*self.data.omega(realtime)**2) * 
                        self.data.leftposL() * self.data.rho_f() + fluxmin[0]
                    )
        #print 'sdot', u_t[0]

        return u_t

    def save_output(self, dir='test'):
        """save for every outputtime a file with the data
        """
        if not os.path.isdir(dir):
            print "Error, save not possible, directory",dir," does not exist"
            sys.exit()
            
        i = 0
        waterheight = zeros( len(self.timeval), float)
        timeout =  zeros( len(self.timeval), float)
        for time in self.timeval : 
            if mod(int(100*time),100*self.data.outputtime()) == 0:
                filename=os.path.join(dir,dir+'_solat_'+(time+self.timeinit).__str__())
                outfile=open(filename,'w')
                outfile.write('# %g %i  # timesec nrlines\n' % 
                        (time+self.timeinit, len(self.grid_x)))
                outfile.write('#pos    pressure\n')
                waterheight[i] = self.sol[i][0]
                timeout[i] = time+self.timeinit
                for (grid, val) in zip(self.grid_x,self.sol[i][1:]) :
                    outfile.write('%g %g \n' %(grid, val))
                outfile.close()
            i=i+1
        sampleLeft = self.data.leftposL()
        containerLeft = self.data.leftpos()
        filename=os.path.join(dir,dir+'_watertable')
        outfile=open(filename,'w')
        outfile.write('#time  height  pressure_left\n')
        for time, pos in zip(timeout, waterheight):
            if not pos ==  0.:
                pres = self.calc_presL(pos, time)
                outfile.write('%g %g %g \n' %( 
                        time, sampleLeft - (containerLeft +pos), pres) )
        outfile.close()
