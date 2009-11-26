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

from LinDiff import DiffusionPDE
import numpy as N
import scipy as S
import UtilsBm
from scipy import zeros, ones, array,  arange,  mod
from scipy.integrate import dae,  odeint
import copy, os

#-------------------
#
# module constant
#
#-------------------

g = 9.80665 #standard gravity, m/s^2
rho_f = 998.2071 # density at 20 degree Celcius, kg / m^3

#-------------------
#
# Solver class
#
#-------------------

class SatUnsatFlow(DiffusionPDE):
    """Model saturated - unsaturated flow in a sample with possible
        rotational force
    input:  grid: the grid between [0,1]
                uvalinit : initial values over the grid
                data : Data of the sample
        Solution method: ....
    """
    def __init__(self, grid, uvalinit, data, dir='test'):
        self.data = data
        DiffusionPDE.__init__(self, dir, None)

        self.grid = grid.copy()
        self.lx = len(self.grid)
        self.grid_cells = zeros(self.lx-1,  float)
        self.grid_cells = (self.grid[1:] + self.grid[:-1])/2.
        
        #fixed grid data
        self.alpha = zeros(self.lx-1, float)
        self.alpha[:] = self.data.L()*(self.grid[1:] - self.grid[:-1])
        
        #time data
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
        
        #solving method
        self.solmeth = 'finitevolume_ode'
        self.lsodaprecr = 1e-5
        self.lsodapreca = 1e-4
        
        #calc pos change in porous stone to sample
        self.leftpos_por1 = self.data.leftposL()
        self.leftpos_samp = self.data.leftposL() + self.data.th_p1()
        self.leftpos_por2 = self.data.leftposL() + self.data.L() \
                                - self.data.th_p2()

        #returning parameters
        self.eps = self.data.eps()
        self.eps_th_p = self.eps / (self.data.theta_s_p() - self.data.theta_r_p())
        self.eps_th = self.eps / (self.data.theta_s() - self.data.theta_r())
        self.conductivity_sat_p = self.data.k_p() * rho_f * g / self.data.mu()
        self.conductivity_sat_s = self.data.k() * rho_f * g / self.data.mu()
        #van genuchten
        self.n_s = self.data.ngen()
        self.n_p = self.data.ngen_p()
        self.m_s = 1.-1./self.n_s
        self.m_p = 1.-1./self.n_p

        self.uvalinit = copy.deepcopy(uvalinit)
        ##input is in pressure, convert to head:
        ##self.uvalinit = zeros(len(uvalinit), float)
        ##print uvalinit;
        ##for i in range(len(uvalinit)):
        ##    if uvalinit[i] >= 0.:
        ##        if self.data.leftposL() + self.data.L() * self.grid[i] < self.leftpos_samp:
        ##            self.uvalinit[i] = 1. + self.eps_th_p * uvalinit[i]/(rho_f*g)
        ##        elif self.data.leftposL() + self.data.L() * self.grid[i] < self.leftpos_por2:
        ##            self.uvalinit[i] = 1. + self.eps_th * uvalinit[i]/(rho_f*g)
        ##        else:
        ##            self.uvalinit[i] = 1. + self.eps_th_p * uvalinit[i]/(rho_f*g)
        ##    else:
        ##        if self.data.leftposL() + self.data.L() * self.grid[i] < self.leftpos_samp:
        ##            self.uvalinit[i] = self.saturation_h(uvalinit[i]/(rho_f*g), 0)
        ##        elif self.data.leftposL() + self.data.L() * self.grid[i] < self.leftpos_por2:
        ##            self.uvalinit[i] =self.saturation_h(uvalinit[i]/(rho_f*g))
        ##        else:
        ##            self.uvalinit[i] = self.saturation_h(uvalinit[i]/(rho_f*g), 2)
        #data storage
        self.head = zeros(self.lx,  float)
        self.sat = zeros(self.lx,  float)
        self.conductivity = zeros(self.lx,  float) # in the base grid (=edges)
        self.flux_left_head = zeros(self.lx,  float) # flux left of cell on edge
        #calculation of some parameters we will need:
        self.par_sat_eq = zeros(self.lx-1,  float)
        for i in range(self.lx-1):
            if self.data.leftposL() + self.data.L() * self.grid_cells[i] < self.leftpos_por1:
                self.par_sat_eq[i] = 1./(self.data.theta_s_p() - 
                                                self.data.theta_r_p())
            elif self.data.leftposL() + self.data.L() * self.grid_cells[i] <= self.leftpos_samp:
                self.par_sat_eq[i] = 1./(self.data.theta_s() - 
                                                self.data.theta_r())
            else:
                self.par_sat_eq[i] = 1./(self.data.theta_s_p() - 
                                                self.data.theta_r_p())
        self.testje = 0

    def calc_headL(self, s = 0., t=0.):
        """
        Calculate the head to the left given s, the position of the water
         table relative to begin of container and t, time of calculation
        We have a contribution of pressure of pure water.
        Head is pressure/(rho_f * g )
        """
        contr = 0.
        rho_om2 = (self.data.omega(t) ** 2) / 2. /g 
        pos_s = self.data.leftpos() + s
        posL = self.data.leftposL()

        #we do not include under pressure!
        #assert (pos_s < posL),  'pos_s : %g, left sample : %g' % (pos_s, posL)
        if (pos_s > posL):
            print  'WARNING, no water: pos_s : %g, left sample : %g' % (pos_s, posL)
        #we have a contribution of pure water
        contr += posL **2 - pos_s **2
        return contr * rho_om2

    def init_cond_fv(self, init):
        """Calculate initial head for finite volume, so in cells
           from init profile in grid which is in head !
        """
        headleft = self.calc_headL(self.data.s0(), 0.)  #C0 of jozef
        headright = self.data.h2()
        #calculate initial heads from given initial cond
        self.head_cells = zeros(self.lx-1,  float)
        self.head_cells[:] = (init[:-1] + init[1:]) / 2.
        #calculate from this initial saturation
        self.sat_cells = zeros(self.lx-1,  float)
        self.lastgridc_por1 = 0
        self.lastgridc_samp = self.lx-2
        self.lastgrid_por1 = 0
        self.lastgrid_samp = self.lx-1
        for i in range(self.lx-1):
            if self.head_cells[i] >= 0.:
                if self.data.leftposL() + self.data.L() * self.grid_cells[i] < self.leftpos_samp:
                    self.sat_cells[i] = 1. + self.eps_th_p * self.head_cells[i]
                elif self.data.leftposL() + self.data.L() * self.grid_cells[i] < self.leftpos_por2:
                    self.sat_cells[i] = 1. + self.eps_th * self.head_cells[i]
                else:
                    self.sat_cells[i] = 1. + self.eps_th_p * self.head_cells[i]
            else:
                if self.data.leftposL() + self.data.L() * self.grid_cells[i] < self.leftpos_samp:
                    self.sat_cells[i] = self.saturation_h(self.head_cells[i], 0)
                elif self.data.leftposL() + self.data.L() * self.grid_cells[i] < self.leftpos_por2:
                    self.sat_cells[i] =self.saturation_h(self.head_cells[i])
                else:
                    self.sat_cells[i] = self.saturation_h(self.head_cells[i], 2)
            if self.data.leftposL() + self.data.L() * self.grid_cells[i] < self.leftpos_samp:
                self.lastgridc_por1 = i
            elif self.data.leftposL() + self.data.L() * self.grid_cells[i] <= self.leftpos_por2:
                self.lastgridc_samp = i
            if self.data.leftposL() + self.data.L() * self.grid[i] < self.leftpos_samp:
                self.lastgrid_por1 = i
            elif self.data.leftposL() + self.data.L() * self.grid[i] <= self.leftpos_por2:
                self.lastgrid_samp = i

    def saturation_h(self, head, mat=1):
        """Calculate saturation, if mat=1: the sample values, otherwise the
           values of the porous stone
        """
        if mat == 1: 
            return 1./N.power( 1.+ N.power(N.abs(self.data.alphagen()*head)
                                           , self.n_s)
                              , self.m_s)
        else:
            return 1./N.power( 1.+ N.power(N.abs(self.data.alphagen_p()*head)
                                           , self.n_p)
                              , self.m_p)

    def head_sat(self, sat):
        """Calculate head from a saturation given in cells - vectorized
        """
        saturated = (sat >= 1)
        unsaturated = N.invert(saturated)
        #calculate the saturated depending on wich material
        self.head_cells[:self.lastgridc_por1+1][saturated[:self.lastgridc_por1+1]] = \
            (sat[:self.lastgridc_por1+1][saturated[:self.lastgridc_por1+1]] - 1.)/self.eps_th_p
        self.head_cells[self.lastgridc_samp+1:][saturated[self.lastgridc_samp+1:]] = \
            (sat[self.lastgridc_samp+1:][saturated[self.lastgridc_samp+1:]] - 1.)/self.eps_th_p
        self.head_cells[self.lastgridc_por1+1:self.lastgridc_samp+1][saturated[self.lastgridc_por1+1:self.lastgridc_samp+1]] =\
            (sat[self.lastgridc_por1+1:self.lastgridc_samp+1][saturated[self.lastgridc_por1+1:self.lastgridc_samp+1]] - 1.) / \
                        self.eps_th
        #calculate the unsaturated depending on which material
        inv_np = 1./self.n_p
        inv_mp = 1./self.m_p
        inv_ns = 1./self.n_s
        inv_ms = 1./self.m_s
        self.head_cells[:self.lastgridc_por1+1][unsaturated[:self.lastgridc_por1+1]] = \
                1./self.data.alphagen_p() * \
                    N.power(-1. + 
                        1./N.power(sat[:self.lastgridc_por1+1][unsaturated[:self.lastgridc_por1+1]], 
                                   inv_mp)
                            , inv_np )
        self.head_cells[self.lastgridc_samp+1:][unsaturated[self.lastgridc_samp+1:]] = \
                1./self.data.alphagen_p() * \
                    N.power(-1. + 
                        1./N.power(sat[self.lastgridc_samp+1:][unsaturated[self.lastgridc_samp+1:]], 
                                   inv_mp)
                            , inv_np )
        self.head_cells[self.lastgridc_por1+1:self.lastgridc_samp+1][unsaturated[self.lastgridc_por1+1:self.lastgridc_samp+1]]\
            =   1./self.data.alphagen() * \
                    N.power(-1. + 
                        1./N.power(sat[self.lastgridc_por1+1:self.lastgridc_samp+1][unsaturated[self.lastgridc_por1+1:self.lastgridc_samp+1]], 
                                   inv_ms)
                            , inv_ns )
        if self.testje < 2:
            print 'head_cells',  self.head_cells

    def conductivity_calc(self, sat_edges):
        """Calculate the conductivity in the grid (=edges), given the saturation
            in edges
            Vectorized
        """
        if self.testje <2:
            print 'cond calc 1',  sat_edges
        saturated = (sat_edges >= 1)
        unsaturated = N.invert(saturated)
        #set value in saturated part
        self.conductivity[:self.lastgrid_por1+1]\
                         [saturated[:self.lastgrid_por1+1]] = \
                self.conductivity_sat_p
        self.conductivity[self.lastgrid_samp+1:]\
                         [saturated[self.lastgrid_samp+1:]] = \
                self.conductivity_sat_p
        self.conductivity[self.lastgrid_por1+1:self.lastgrid_samp+1]\
                         [saturated[self.lastgrid_por1+1:self.lastgrid_samp+1]]\
                = self.conductivity_sat_s
        #set value in unsaturted part with van genuchten sat
        self.conductivity[:self.lastgrid_por1+1]\
                         [unsaturated[:self.lastgrid_por1+1]] = \
                self.conductivity_sat_p * N.power(
                        sat_edges[:self.lastgrid_por1+1]\
                                 [unsaturated[:self.lastgrid_por1+1]]
                                                 , 0.5) * \
                    N.power (1- N.power(1- N.power( 
                        sat_edges[:self.lastgrid_por1+1]\
                                 [unsaturated[:self.lastgrid_por1+1]]
                                                   , 1./self.m_p)
                                        , self.m_p)
                            , 2.)
        self.conductivity[self.lastgrid_samp+1:]\
                         [unsaturated[self.lastgrid_samp+1:]] = \
                self.conductivity_sat_p * N.power(
                        sat_edges[self.lastgrid_samp+1:]\
                                 [unsaturated[self.lastgrid_samp+1:]]
                                                 , 0.5) * \
                    N.power (1- N.power(1- N.power( 
                        sat_edges[self.lastgrid_samp+1:]\
                                 [unsaturated[self.lastgrid_samp+1:]]
                                                   , 1./self.m_p)
                                        , self.m_p)
                            , 2.)
        self.conductivity[self.lastgrid_por1+1:self.lastgrid_samp+1]\
                         [unsaturated[self.lastgrid_por1+1:self.lastgrid_samp+1]] = \
                self.conductivity_sat_p * N.power(
                        sat_edges[self.lastgrid_por1+1:self.lastgrid_samp+1]\
                                 [unsaturated[self.lastgrid_por1+1:self.lastgrid_samp+1]]
                                                 , 0.5) * \
                    N.power (1- N.power(1- N.power( 
                        sat_edges[self.lastgrid_por1+1:self.lastgrid_samp+1]\
                                 [unsaturated[self.lastgrid_por1+1:self.lastgrid_samp+1]]
                                                   , 1./self.m_p)
                                        , self.m_p)
                            , 2.)
        if self.testje <2:
            print 'cond calc 2',  self.conductivity

    def solve(self):
        #bc data must be set
        import time
        e0 = time.time()
        c0 = time.clock()
        self.shifttime = 0.
        #integrate it, jacobian is banded matrix ml=1, mu=1
        # do deepcopy of second par as this is changed internally with new value
        self.sol = self._solve(self.uvalinit, self.timeval)
        #we have a solution with the present parameters
        self.changedparam = False
        elapsed_time = time.time()- e0
        cpu_time = time.clock() - c0
        print ' Solved', self.__class__.__name__, ' in', \
                elapsed_time, 'seconds, cpu time', \
                cpu_time
        print self.sol

    def _solve(self, init, timeval):
        """Do necessary to interface different solvers"""
        if self.solmeth == 'finitevolume_ode' :
            init_sat = self.init_cond_fv(init)
            initode = zeros(self.lx, float)
            initode[0] = self.data.s0()
            initode[1:] = (copy.deepcopy(self.sat_cells))[:]
            print 'init',  initode
            sol = odeint(self.f_cons_fv, initode, timeval,
                            rtol=self.lsodaprecr, atol=self.lsodapreca, 
                            printmessg = 1, 
                            ml=1, mu=1,  h0=1e-7)
            return sol
        else :
            print 'ERROR : unknown solution method to solve ODE system'
            sys.exit()
    
    def f_cons_fv(self,u,t) :
        ''' 
            right handed side of 1st - Order ODE system modelling the 
            saturated unsaturated flow, u[0]=s, u[1:]=saturation
                u_t = 1/(th_s-th_r) * diff_r (K * h_r - omega^2 * r/g)
            with u saturation if h<0, and if h>0 the artificial: 
                        u=1+eps/(th_s-th_r) h
            We work with cell data, all cells are full cells
        '''
        # u is cell data in self.grid_cells, at time t we need to return u_t
        u_t    = zeros(self.lx, float)
        #clean up input, negative saturation is not allowed! 
        if self.testje <4:
            print 'u', u
        tmp = u<1e-12
        u[1:][tmp[1:]] = 1e-12
        # first, from u[0], compute head at left border:
        self.head[0] = self.calc_headL(u[0], t)
        # compute head from saturation in experiment
        self.head_sat(u[1:])
        #compute head in edges, so over the starting grid
        self.head[1:-1] = (self.head_cells[:-1] + self.head_cells[1:])/2.
        self.head[-1] = self.data.h2()  #dirichlet for now, should be improved
        #compute saturation in edges
        self.sat[1:-1] = (u[1:-1] + u[2:])/2.
        if self.data.th_p1() > 0.:
            if self.head[0] >= 0.:
                if self.testje == 0:
                    print 'yes',  self.eps_th_p * self.head[0]
                self.sat[0] = 1. + self.eps_th_p * self.head[0]
            else:
                self.sat[0] = self.saturation_h(self.head[0], 0)
        else:
            if self.head[0] >= 0.:
                self.sat[0] = 1. + self.eps_th * self.head[0]
            else:
                self.sat[0] = self.saturation_h(self.head[0], 1)
        if self.data.th_p2() > 0.:
            if self.head[-1] >= 0.:
                self.sat[-1] = 1. + self.eps_th_p * self.head[-1]
            else:
                self.sat[-1] = self.saturation_h(self.head[-1], 0)
        else:
            if self.head[-1] >= 0.:
                self.sat[-1] = 1. + self.eps_th * self.head[-1]
            else:
                self.sat[-1] = self.saturation_h(self.head[-1], 1)
        #compute conductivity in edges
        self.conductivity_calc(self.sat)
        
        om2og = N.power(self.data.omega(t), 2) / g
        #now compute fluxleft and fluxright
        self.flux_left_head[1:-1] = \
                (self.head_cells[1:] -self.head_cells[:-1]) \
                 /((self.alpha[:-1]+self.alpha[1:])/2.)
        ##self.flux_left_head[0] = UtilsBm.deriv131(self.leftpos_por1, self.head[0],
        ##                self.leftpos_por1+self.alpha[0], self.head[1],
        ##                self.leftpos_por1+self.alpha[0]+self.alpha[1], 
        ##                                                        self.head[2])
        self.flux_left_head[0] = (self.head[1]-self.head[0])/self.alpha[0]
        self.flux_left_head[-1] = 0.  #Hom Neumann condition to right
        #change of positon of the watertable is first equation:
        #     d watertab/dt = q (darcy flux)
        #     distance from sensor is first eq, so invert: * -1 !
        if self.testje <2:
            print 'flux left head',  self.flux_left_head
        u_t[0] = -1.*self.conductivity[0] * (self.flux_left_head[0] \
                                            - om2og * self.leftpos_por1)
        #change in saturation is second eqauation to last equation
        u_t[1:] = self.par_sat_eq[:]/ self.alpha[:] * (
                    self.conductivity[1:]*(self.flux_left_head[1:] 
                        - om2og * (self.leftpos_por1 + self.data.L() 
                                                            * self.grid[1:]))
                    - self.conductivity[:-1]*(self.flux_left_head[:-1] 
                        - om2og * (self.leftpos_por1 + self.data.L() 
                                                            * self.grid[:-1]))
                    )
        if self.testje <4:
            print t,  u_t
        self.testje += 1
        return u_t

    def showmass(self,  dir='test'):
        pass

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
                        (time+self.timeinit, len(self.grid_cells)))
                outfile.write('#pos    saturation\n')
                waterheight[i] = self.leftpos_samp - self.sol[i][0] \
                                    -self.leftpos_por1
                timeout[i] = time+self.timeinit
                for (grid, val) in zip(self.grid_cells,self.sol[i][1:]) :
                    outfile.write('%g %g \n' %(grid, val))
                outfile.close()
            i=i+1
        filename=os.path.join(dir,dir+'_watertable')
        outfile=open(filename,'w')
        outfile.write('#time  height \n')
        for time, pos in zip(timeout, waterheight):
            if not pos ==  0.:
                outfile.write('%g %g \n' %( 
                        time, pos) )
        outfile.close()
