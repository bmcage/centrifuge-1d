#!/usr/bin env python

# Copyright (C) 2000-2006  B. Malengier
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
from interpolateBM import interp1dBM
import GridUtils
from diffusion import *

#-------------------------------------------------------------------------
#
# General diffusion PDE class
#
#-------------------------------------------------------------------------

class DiffusionPDE:
    ''' general class of what a diffusion equation should be able to do
    '''
    def __init__(self, outputdir, diffusion):
        self.dir     = outputdir
        self.diff    = diffusion
        self.timeval = array([0.,])
        #on modification of data, changedparam is set to true
        self.changedparam = True
        #standard a method is used for a direct problem
        self.adjoint = False
        self.sourceterm_dt = False
        self.bcvalparam = []
        self.timeinit = 0.
        
    def solve(self):
        ''' solve the equation  '''
        pass
    
    def set_outputdir(self, str) :
        self.dir = str
        
    def get_outputdir(self) :
        return self.dir 
    
    def set_adjoint(self):
        self.adjoint = True
    def set_direct(self):
        self.adjoint = False
    
    def get_diffusion(self):
        return self.diff
    
    def set_diffusion(self,diffusion):
        self.diff = diffusion
        
    def timevalues(self):
        '''
        returns a list of all timevalues at which a solution is calculated
        '''
        return self.timeval
    
    def _set_timevalues(self, times):
        '''
        set a list of all timevalues at which a solution is calculated
        '''
        self.timeval = times
 
    def set_boundary_data(self, bcL, bcR, bcvalL=0., bcvalR=0., 
                            bcvalident=True, comp=1):
        """ Set easy normal boundary flux data for 1 component
            Why flux? We solve with finite volume, so flux over edge needed
        """
        #set boundary conditions
        if comp == 1:
            self.bcL_u = bcL
            self.bcR_u = bcR
            indexstart = 0
            self.bcvalparam = []
        else:
            self.bcL_v = bcL
            self.bcR_v = bcR
            indexstart = len(self.bcvalparam)
        
        self.bcvalparam.append(bcvalL)
        if bcL =='homNeu' : self.bcvalparam[-1] = 0.
        if bcL =='homNeu' or bcL =='Neu' :
            #given time and concentration on border return du/dx
            if comp == 1:
                self.bcvalL_u = lambda t,u : self.bcvalparam[indexstart] #function given time,conc return Neumann cond left
            else:
                self.bcvalL_v = lambda t,u : self.bcvalparam[indexstart] #function given time,conc return Neumann cond left
        elif bcL == 'homRob' :
            #given time and concentration on border return du/dx
            if comp == 1:
                self.bcvalL_u = lambda t,u : self.bcvalparam[indexstart] * u
            else:
                self.bcvalL_v = lambda t,u : self.bcvalparam[indexstart] * u
        else : print 'ERROR: Unknown border condition', self.bcL_u;sys.exit()
        if bcvalident :
            indexpar = indexstart
            bcR = bcL
            #overwrite right values with left values
            if comp == 1:
                self.bcR_u = self.bcL_u
            else:
                self.bcR_v = self.bcL_v
        else :
            self.bcvalparam.append(bcvalR)
            indexpar = len(self.bcvalparam)
        if bcR =='homNeu' : self.bcvalparam[indexpar] = 0.
        if bcR =='homNeu' or bcR =='Neu' :
            #given time and concentration on border return du/dx
            if comp == 1:
                self.bcvalR_u = lambda t,u : self.bcvalparam[indexpar] #function given time,conc return Neumann cond left
            else:
                self.bcvalR_v = lambda t,u : self.bcvalparam[indexpar] #function given time,conc return Neumann cond left
        elif self.bcR_u == 'homRob' :
            #given time and concentration on border return du/dx
            if comp == 1:
                self.bcvalR_u = lambda t,u : self.bcvalparam[indexpar] * u
            else:
                self.bcvalR_v = lambda t,u : self.bcvalparam[indexpar] * u
        else : print 'ERROR: Unknown border condition', self.bcR_u;sys.exit()

    def set_boundary_data_MC(self, bcL, bcR, bcL_v, bcR_v, bcvalL=0., bcvalR=0., 
                            bcvalL_v=0., bcvalR_v=0., bcvalident=True):
        """ Set easy normal boundary flux data for 2 components
            Why flux? We solve with finite volume, so flux over edge needed
        """
        #set boundary conditions first comp
        self.set_boundary_data(bcL, bcR, bcvalL, bcvalR, bcvalident, comp=1)
        #set boundary conf second comp
        self.set_boundary_data(bcL_v, bcR_v, bcvalL_v, bcvalR_v, bcvalident, 
                               comp=2)

    def set_boundary_data_func(self, func, pos=0, comp=1):
        """ Set boundary flux data for 1 component where factor is a function
            in time and concentration. 
            So this can be used for eg  partial_x u = function(t) * u, in 
                which case func(t,u) = function(t) * u
            If pos=0 then left border, otherwise right
            Why flux? We solve with finite volume, so flux over edge needed
        """
        if pos == 0:
            if comp == 1:
                self.bcvalL_u = lambda t,u : func(t, u)
            else:
                self.bcvalL_v = lambda t,u : func(t, u)
        else:
            if comp == 1:
                self.bcvalR_u = lambda t,u : func(t, u)
            else:
                self.bcvalR_v = lambda t,u : func(t, u)
        #make sure bc data is recognized as set
        self.bcvalparam = []

    def set_time_data(self,endtime, timestep, nrsteps, type='timestep', 
                              time_experiments=[],expmayadd=True):
        ''' We solve in interval [0,endtime]. 
            We need to keep output data: 
               if type= timestep: every timestep seconds
               if type= nrsteps: nrsteps intervals.
           time_experiments are the times at which experiments happen
           if expmayadd, then adding these is allowed, 
              if expmayadd == False, then the experimental timesteps should be in the 
                              default times. (needed if equidistant timesteps are needed)
        '''
        self.end = endtime
        if type == 'timestep' :
            assert(timestep> 0.)
            self.timestep = timestep
            nrsteps_fl  = self.end / float(timestep)
            nrsteps_int = int(nrsteps_fl)
            #we need to arrive at end with integer number of steps:
            assert(nrsteps_fl-nrsteps_int == 0.)
            self.nrsteps = nrsteps_int
        elif type == 'nrsteps' :
            assert(int(nrsteps))
            self.nrsteps = int(nrsteps)
            self.timestep = self.end / self.nrsteps
        
        timeval = []
        for i in arange(self.nrsteps+1):
            timeval.append(i*self.timestep)
        
        lenorig = len(timeval)
        if not(time_experiments == []) :
            timeval = sorted(UtilsBm.merge_nodup(timeval, time_experiments))
            if not expmayadd :
                #when storing, we need that experiments are during recorded timesteps
                assert(len(timeval)==lenorig)
        self._set_timevalues(array(timeval))
        
    def set_source_dt(self, timeval, fval):
        '''source term is a delta(t-t_i)f_i(x) function, in which 
            t_i = timeval[i] and f_i(x_j) = fval[i][j] with x_j the grid of the SOLUTION !!
        '''
        assert(len(timeval) == len(fval))
        
        if not timeval == [] :
            self.sourceterm_dt = True
            #from self.timeval we make pieces
            self.time_source_dt = []
            sourcetime = timeval[0]
            ind0 = 0
            time0=0.
            for time in timeval :
                ind = list(self.timeval).index(time)
                self.time_source_dt += [[x-time0 for x in list(self.timeval[ind0:ind+1])]]
                ind0 = ind+1
                time0 = time
            self.time_source_dt += [[x-time0 for x in list(self.timeval[ind0:])]]
            self.fval_source_dt =fval
        else :
            self.sourceterm_dt = False
    
    def plot_sol(self, dir='test'):
        '''
        plot the solution on screen. This will depend on how solve returns and stores
        the solution
        '''
        pass
        
    def mass(self):
        """
        calculate the mass from the concentrations in the solution. 
        should return a list of data over all timevalues
        """
        pass
    
    def mass_calcmethod(self):
        '''
        return how mass is calculated
        '''
        pass
    
    def showmass(self, dir='test', unit='[10^-3 mol/mm^2]') :
        """
        print the mass balance nicely out on screen
        """
        try:
            mass = self.mass()
            timeval = self.timevalues()
            g   = Gnuplot.Gnuplot(persist=0)
            p1=Gnuplot.Data(timeval, mass, with_='lines')
            g.title('Mass balance ('+ self.mass_calcmethod()+')')
            g.xlabel('t [s]')
            g.ylabel('line Mass '+unit)
            g.plot(p1)
            if not os.path.isdir(dir):
                print ("Error, plot not possible, directory",dir," does not exist")
                sys.exit()
            g.hardcopy(filename=os.path.join(dir,dir+'_massbal.ps'), enhanced=1,mode='eps',
                    color=1,fontname='Times-Roman',fontsize=28)
        except : 
            print 'ERROR : showmass() showing mass balance failed'

    def residuals(self) :
        '''Compute the difference between experiment and solution,
            over deviation
            return as an array
        '''    
        return None  #not implemented

class Lin1DDiffusion(DiffusionPDE):
    '''implement 1 D linear diffusion u_t = (D(x) u_x)_x 
        input:  grid: the grid between [0,1]
                length : real length of the grid
                left_cutoff : how much difference between 0 and start grid
                uvalinit : initial values over the grid
                diff  : a diffusion object (D(x) returns diffusion)
        Solution method: ....
    '''
    def __init__(self,grid,uvalinit,diff,dir,length_grid=1.,left_cutoff=0.,
                    solmeth='conservative_fv'):
        DiffusionPDE.__init__(self, dir, diff )
        self.grid = grid
        #np   = len(grid)
        self.length = length_grid
        self.cutoff = left_cutoff
        self.uvalinit = uvalinit
        
        self.sol = None
        self.bcvalparam = None
        
        # new object, so changed parameters
        self.changedparam = True
        
        self.solmeth=solmeth
        
        # set data needed to solve the problem
        #set lsoda band matrix and precision
        self.lsodaprec = 1e-14
        self.lsodajtype= 5
        #needed varibles, depending on solution method, only calcuate once
        if self.solmeth == 'conservative_fv' :
            #finite volumes: cell values are the init cond
            #uvalinit_fv = zeros(np-1,float)
            self.uvalinit_fv = (uvalinit[:-1] + uvalinit[1:])/2.
            #self.Deltamin = zeros(np-1,float)
            self.Deltamin = self.grid[1:]-self.grid[:-1]
            #real grid where the solution exists on
            #solution is in the midpoints
            self.grid_x = (self.grid[:-1] + self.grid[1:]) / 2.
            self.grid_x = self.cutoff + self.grid_x *self.length
            #we need the real grid too every time we calculate
            self.grid_actually = self.cutoff + self.grid *self.length
        else :
            print 'solution method not implemented, LinDiff.py'
            sys.exit()

    def set_init_data(self,uvalinit,celldata=False):
        ''' initial data over the grid.
            Note: for FV approx, this will be converted to cell data!
            if celldata=True, the uval is over cells, not in grid
        '''
        if celldata :
            #we should not be given celldata if not fin vol method
            #as we do not know how to extrapolate to boundary!
            assert(self.solmeth == 'conservative_fv')
            self.uvalinit_fv = uvalinit
        else :
            self.uvalinit = uvalinit
            if self.solmeth == 'conservative_fv' :
                #finite volumes: cell values are the init cond
                self.uvalinit_fv = (self.uvalinit[:-1] + self.uvalinit[1:])/2.
        self.changedparam = True
        
    def _solve(self, init, timeval):
        if self.solmeth == 'conservative_fv' :
             return odeint(self.f_cons_fv, copy.deepcopy(init),timeval, \
                            rtol=self.lsodaprec, atol=self.lsodaprec, printmessg = 1, \
                            ml=1, mu=1)
        else :
            print 'ERROR : unknown solution method to solve ODE system'
            sys.exit()
        
    def solve(self):
        #bc data must be set
        assert(self.bcvalparam is not None)
        import time
        e0 = time.time()
        c0 = time.clock()
        np   = len(self.grid)
        self.shifttime = 0.
        #integrate it, jacobian is banded matrix ml=1, mu=1
        # do deepcopy of second par as this is changed internally with new value
        if self.sourceterm_dt :
            #split time interval in pieces, run solution in piece, and add 
            #  sourceterm before starting next step.
            i = 0
            init = self.uvalinit_fv
            self.sol = zeros((len(self.timeval),len(self.grid_x)),float)
            prevlen = 0
            l = len(self.time_source_dt)-1 
            for times in self.time_source_dt :
                self.sol[prevlen:prevlen+len(times)] = self._solve(init, times)
                # set new init data
                if i < l :
                    init = array(self.sol[prevlen+len(times)-1][:])
                    init[:] = init[:] + array(self.fval_source_dt[i])[:]
                else :
                    init = None
                prevlen = prevlen+len(times)
                #internally in self._solve we need to know the REAL time!
                self.shifttime += times[-1]
                i = i+1
        else :
            self.sol = self._solve(self.uvalinit_fv, self.timeval)
        #we have a solution with the present parameters
        self.changedparam = False
        elapsed_time = time.time()- e0
        cpu_time = time.clock() - c0
        print ' Solved Lin1DDiffusion in', elapsed_time, 'seconds, cpu time', \
                cpu_time
    
    def f_cons_fv(self,u,t) :
        ''' right handed side of 1st - Order ODE system 
                conservative form: u_t = flux_right - flux_left
                                       = 1/Delta * (D(x)u_x)_+ - (D(u)u_x)_-
           We work with cell data, all cells are full cells
        '''
        np   = len(self.grid)
        #it happens that the odeint is stopped and restarted, we need to have correct time:
        realtime = self.shifttime + t
        #set all inner points, uvalinit also np-1 large!
        #u_t    = zeros(np-1, float)
        # an adjoint problem internally works with an inverted time
        # calculate the real time.
        #time = t
        #if self.adjoint :
        #    time = self.end - t
            
        #diffusion in cells, needed diffusion on the edges = in gridpoints
        diff   = zeros(np, float)
        diff[0] = 1.
        #only determine diffusion in the internal points (grid_x is internal!)
        diff[1:-1]   = self.diff.D(realtime, self.grid_actually[1:-1])
        diff[-1] = 1.
        #fluxes on all edges (=correspond to gridpoints)
        fluxmin = zeros(np,float)
        fluxmin[1:-1] = (u[1:]-u[:-1])/     \
                        ((self.Deltamin[:-1]+self.Deltamin[1:])/2.)
        #border flux is the border Neumann cond
        uLborder = u[0] #first approx, value border is cell value
        fluxleftBC = self.bcvalL_u(realtime,uLborder) 
        #normal to the left is -1!
        fluxmin[0] = -1*fluxleftBC
        uRborder = u[-1] #first approx, value border is cell value
        fluxrightBC = self.bcvalR_u(realtime,uRborder)
        fluxmin[-1] = fluxrightBC
        
        t1     = self.length * self.length
        t2     = 1./t1
        # u_t = (Flux to right - Flux to left) / size_cell
        u_t = t2 * (diff[1:]*fluxmin[1:]- \
                            diff[:-1]*fluxmin[:-1]) \
                        / (self.Deltamin[:])
        return u_t       

    def solution(self, time):
        '''return the solution at a required time, use only for adjoint method as
            only then we require nice equidistant timeslices, which this method uses !! '''
        TODO # this method needs to be tested and updated as in NonLinDiff.py
        post =  int(time / float(self.nrsteps))
        if post == len(self.timeval) :
            post = post -1 
        #small test:
        step = self.timeval[1]-self.timeval[0]
        if not (post * step <= time and (post+1)*step >= time) :
            print 'error in def solution NonLinDiff.py '
            sys.exit
        #note that this solution is NOT over grid, but over grid_x, used in the solver!
        return self.sol[post+1] * (time-post*step) / step + \
                 self.sol[post] * (1.- (time-post*step)/step)
