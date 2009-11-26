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
import Gnuplot
import pylab as pl
import copy

#-------------------------------------------------------------------------
#
# local modules
#
#-------------------------------------------------------------------------

import UtilsBm
from UtilsBm import piecewise3
from interpolateBM import interp1dBM
import GridUtils

#-------------------------------------------------------------------------
#
# A constant diffusion 
#
#-------------------------------------------------------------------------
 
class Diffusion_const:
    ''' A constant diffusion
    '''
    def __init__(self, value):
        self.param = [value]
        self.nrparam = 1    
        self.type = 'Function'
        self.modified = True
        
        self.make_diff_func()

    def make_diff_func(self) :
        ''' everytime parameters are changed, the diffusion 
            function must be REMADE !! 
        '''
        self.diff  = lambda y: zeros(len(y), float) + self.param[0]
        self.diff_u= lambda y: zeros(len(y), float)
        #value of diffusion is again in agreement with par
        self.modified = False

    def nrParam(self) :
        return self.nrparam
    
    def get_param(self, number) :
        return self.param[number]
    
    def set_param(self, number, val) :
        self.param[number] = val
        self.modified = True
        
    def D(self,v):
        ''' function value of diffusion in array v'''
        if self.modified : 
            self.make_diff_func()
            print ' test par dif ' , self.param
        return self.diff(v)
    
    def D_u(self,v):
        ''' function value of derivative of diffusion to u in array v'''
        if self.modified : self.make_diff_func()
        if self.diff_u :
            return self.diff_u(v)
        else:
            return None

#-------------------------------------------------------------------------
#
# A diffusion depending on x via a grid function u: D(x) = D(u(x))
#
#-------------------------------------------------------------------------
 
class Diffusion_xt_uxt:
    ''' a diffusion function, returning given the time t and position x, 
        the diffusion D(x,t), using a function u, so internally D(u) can 
        be calculated, knowing u(x,t) gives D(u(x,t)) 
    '''
    def __init__(self, gridx, deltax, griduval = True, equidistant=True):
        '''
        set the diffusion over a grid x, which is a numpy list.
        The list should be equidistant, then deltax is used
        If griduval = True : data of u is in the gridpoints
                     = False: data of u is in the cell midpoint
        '''
        self.grid = gridx
        self.lengrid = len(gridx)
        self.stepgrid = deltax
        self.diff = None
        self.griduval = griduval
        self.uval_time = None
        #only equidistant implemented
        assert(equidistant)
        #only cell data implemented
        assert(not(griduval))
        
        
    def set_diff_func(self, diff):
        '''set the diffusion function in u used, eg a Diffusion_u class object
        '''
        self.diff = diff
        
    def set_grid_func_time_method(self,method):
        ''' set the method with wich uxt at a required time can be obtained,
            so method(time) returns an uval
        '''
        self.uval_time_method = method
        
    def set_grid_func(self, uval, time, griduval=True):
        ''' diffusion will be returned based on the values uval, which are values
            correct for time
            If griduval, the given uval is in gridpoints, if False it is in cell
        '''
        if self.griduval == griduval:
            assert(len(uval) == len(self.grid))
        else :
            assert(len(uval) == len(self.grid)-1)
        self.uval = uval
        self.uval_time = time
        if self.griduval == griduval:
            #the uval values corresponding to self.grid
            self.uval_grid = uval
        else :
            print 'diffusion.py-set_grid_func, not implemented '
            sys.exit()
        
    def D(self,time, xval):
        ''' return value of diffusion is (numpy) array of points xval
        '''
        #determine the x coord needed, following are normal slow arrays as result
        if not(time == self.uval_time) :
            #we need to get uval at time 
            self.set_grid_func(self.uval_time_method(time), time, self.griduval)
            
        def indexcorrect(index) :
            i=int(index)
            if i == self.lengrid-1 :
                i -= 1
            return i
        
        
        #we use the fact the grid is equidistant
        xpos = map(indexcorrect,(xval - self.grid[0])/self.stepgrid)
        uval = array(map(lambda x: self.uval_grid[x], xpos))
        uval_right = array(map(lambda x: self.uval_grid[x+1], xpos))
        grid = array(map(lambda x: self.grid[x], xpos))
        grid_right = array(map(lambda x: self.grid[x+1], xpos))
        #return interpolated value of D in x, so u is with (x_i,u_i) and  (x_i+1,u_i+1)
        uvalgrid = (uval_right -uval)/(grid_right-grid)*(xval-grid) + uval

        return self.diff.D(uvalgrid)

#-------------------------------------------------------------------------
#
# A general diffusion depending on u, D(u)
#
#-------------------------------------------------------------------------
    
class Diffusion_u:
    ''' a diffusion function, returning given the value u, the diffusion
        D(u), u being a function from the position x (eg the concentration)
    '''
    def __init__(self,options):
        #self.diffparam = options.diffparam 
        self.splinedata = None
        self.difftype = options.diffusion
        self.modified = True
        #map the given name for diffusion to a function
        if options.diffusion == 'const':
            valhulp = options.diffparam[0].replace('[',' ').replace(']',' ').split(',')
            self.param = [float(x) for x in valhulp]
            self.nrparam = 1    
            self.type = 'Function'
        elif options.diffusion == 'power':
            valhulp = options.diffparam[0].replace('[',' ').replace(']',' ').split(',')
            self.param = [float(x) for x in valhulp]
            self.nrparam = 4
            self.type = 'Function'
        elif (options.diffusion == 'bspline' or 
                options.diffusion == 'bspline1storder' or 
                options.diffusion =='1D2pint'):
            t1 = [[],[]]
            t1[0] =options.diffparam[0].replace('[',' ').replace(']',' ').split(',')
            t1[1] =options.diffparam[1].replace('[',' ').replace(']',' ').split(',')
            self.paramuval = [float(x) for x in t1[0]]
            self.param = [float(x) for x in t1[1]]
            self.nrparam = len(self.param)
            if options.diffusion =='1D2pint' :
                self.type = 'Local'
            else :
                #a bspline is like a localized wavelet, but not completely
                self.type = 'Function'
        elif options.diffusion == 'function' :
            self.param = []
            self.nrparam = len(self.param)
            self.diff = eval(options.diffparam[0])
            self.type = 'Function'
        else:
            print ('Unknown diffusion type given',options.diffusion)
            sys.exit()
        try :
            self.make_diff_func()
        except IndexError:
            print ('Too few parameters given for diffusion type', options.diffusion)
            sys.exit()
            
    def make_diff_func(self) :
        ''' everytime parameters are changed, the diffusion 
            function must be REMADE !! 
        '''
        if self.difftype == 'const':
            self.diff  = lambda y: zeros(len(y), float) + self.param[0]
            self.diff_u= lambda y: zeros(len(y), float) + 0.             
        elif self.difftype == 'power':
            self.diff  = lambda y: y**self.param[0] \
                    * (self.param[1]+self.param[2]*y+self.param[3]*y**2) 
            self.diff_u= lambda y: self.param[0]*y**(self.param[0]-1.) \
                    * (self.param[1]+self.param[2]*y+self.param[3]*y**2) + \
                    y**self.param[0] * (self.param[2]+2.*self.param[3]*y)
        elif self.difftype == 'bspline':
            #create an interp object
            from scipy.interpolate.fitpack import splrep,splev
            #paramuval should be list of u values
            #param should be list of D(u) values
            #create the data of the cubic bspline:
            # no smoother so s = 0
            self.splinedata = splrep(self.paramuval, self.param,
                    xb = None, xe = None, s=0,
                    k = 3, full_output = 0, quiet = 1)
            self.diff   = lambda y: splev(y, self.splinedata, der = 0)
            self.diff_u = lambda y: splev(y, self.splinedata, der = 1)
        elif self.difftype == 'bspline1storder':
            #create an interp object
            from scipy.interpolate.fitpack import splrep,splev
            #paramuval should be list of u values
            #param should be list of D(u) values
            #create the data of the linear bspline:
            # no smoother so s = 0
            self.splinedata = splrep(self.paramuval, self.param,
                    xb = None, xe = None, s=0,
                    k = 1, full_output = 0, quiet = 1)
            self.diff   = lambda y: splev(y, self.splinedata, der = 0)
            self.diff_u = lambda y: splev(y, self.splinedata, der = 1)  
                
        elif self.difftype == '1D2pint':
            #interpolation with 2 points (=piecewise linear)
            self.diff   = GridUtils.GridFunc1D([self.paramuval],self.param)
            self.diff_u = None
        elif self.difftype == 'function' :
            #self.diff is defined when reading the init file
            self.diff_u = None
        else:
            print ('Unknown diffusion type given',options.diffusion)
            sys.exit()
        #value of diffusion is again in agreement with par
        self.modified = False
            
    def nrParam(self) :
        return self.nrparam
    
    def get_param(self, number) :
        return self.param[number]
    
    def set_param(self, number, val) :
        self.param[number] = val
        self.modified = True
        
    def D(self,v):
        ''' function value of diffusion in array v'''
        if self.modified : 
            self.make_diff_func()
            print ' test par dif ' , self.param
        return self.diff(v)
    
    def D_u(self,v):
        ''' function value of derivative of diffusion to u in array v'''
        if self.modified : self.make_diff_func()
        if self.diff_u :
            return self.diff_u(v)
        else:
            return None
            
    def Dsplineminu(self):
        '''return minimum u value in case of spline, None otherwise '''
        if self.splinedata:
            return self.paramuval[0]
        return None
    def Dsplinemaxu(self):
        '''return minimum u value in case of spline, None otherwise '''
        if self.splinedata:
            return self.paramuval[-1]
        return None
    
    def map(self, uval):
        ''' In the case of localised wavelets, we can map a uval to the values 
            in the given wavelet nodes u_ij
        '''
        assert(self.type == 'Local')
        return self.diff.map(uval)
    
    def plot(self,uL, uR, step=0.01, show=False, dir='test'):
        ''' plot D and diff(D,u) between uL and uR 
            and save the plot to a .png file. The dir must exist
        '''
        pl.clf()
        if not os.path.isdir(dir):
            print ("Error, plot diffusion not possible, directory does not exist")
            sys.exit()
        u = pl.arange(uL, uR+step/2., step)
        Du = self.D(u)
        print 'plot D', Du, type(Du)
        pl.plot(u, Du)
        pl.xlabel('Concentration [mol/mm^3]')
        pl.ylabel('Diffusion [1E-12 m^2/s]')
        #pl.title('Diffusion')
        pl.grid(True)
        pl.savefig(os.path.join(dir,dir+'diff'))
        if show:
            pl.show()
        Du_u = self.D_u(u)
        if Du_u == None :
            return
        pl.plot(u, Du_u)
        pl.xlabel('Concentration [mol/mm^3]')
        pl.ylabel('First derivative of Diffusion')
        #pl.title('Diffusion')
        pl.grid(True)
        pl.savefig(os.path.join(dir,dir+'diff_u'))
        if show:
            pl.show()


#-------------------------------------------------------------------------
#
# A general diffusion depending on u and v, D(u, v)
#
#-------------------------------------------------------------------------
    
class Diffusion_u_MC:
    ''' a multicomponent diffusion function, 
        returning given the value u,v the diffusion
        D(u,v), u and v being functions from the position x (eg the concentration)
    '''
    def __init__(self,options):
        self.modified = True
        self.pargrid = False
        self.param = [[],[],[],[],[],[]]  # u_i, v_i, D(u_i,v_i) 
        self.splinedata = None
        self.arrayaware = True
        
        self.difftype = options.diffusion
        self.diff11 = None;self.diff12 = None
        self.diff21 = None;self.diff22 = None;
        self.diff11_u = None;self.diff12_u = None
        self.diff21_u = None;self.diff22_u = None;
        self.diff11_v = None;self.diff12_v = None
        self.diff21_v = None;self.diff22_v = None;
        #map the given name for diffusion to a function
        try:
            if options.diffusion == 'function':
                self.arrayaware = False
                self.param = []
                self.nrparam = len(self.param)
                self.diff11 = eval(options.diffparam[0])
                self.diff12 = eval(options.diffparam[1])
                self.diff21 = eval(options.diffparam[2])
                self.diff22 = eval(options.diffparam[3])
                temp = eval(options.diffparam[4])
                self.umin = temp[0]; self.umax = temp[1]
                self.vmin = temp[2]; self.vmax = temp[3]
            if options.diffusion.split('_')[-1]=='ext1D' :
                (uval, tD11, tD12, vval, tD21, tD22, nru, nrv) = \
                        self._extract_par(options)
                # we extend D11 and D12 value in the v direction
                # we extend D21 and D22 value in the u direction
                for i in range(0,nru) :
                    for j in range(0,nrv) :
                        self.param[0] = self.param[0] + \
                            [uval[i]]
                        self.param[1] = self.param[1] + \
                            [vval[j]]
                        self.param[2] = self.param[2] + \
                            [tD11[i]]
                        self.param[3] = self.param[3] + \
                            [tD12[i]]
                        self.param[4] = self.param[4] + \
                            [tD21[j]]
                        self.param[5] = self.param[5] + \
                            [tD22[j]]   
                self.nrparam = len(self.param[2]) 
                
            if options.diffusion == 'MC_interp2d_ext1D':
                self.param[2] = []
                self.param[3] = []
                self.param[4] = []
                self.param[5] = []
                self.nrparam = 0
                for i in range(0,nru) :
                    td11 = []; td12 = []; td21 = []; td22 = []
                    for j in range(0,nrv) :
                        td11 = td11 + [tD11[i]]
                        td12 = td12 + [tD12[i]]
                        td21 = td21 + [tD21[j]]
                        td22 = td22 + [tD22[j]]   
                    self.param[2].append(td11)
                    self.nrparam += len(td11)
                    self.param[3].append(td12)
                    self.param[4].append(td21)
                    self.param[5].append(td22)
                self.rowlength = len(td11)
                self.pargrid = True
                self.param[0] = uval
                self.param[1] = vval
            elif options.diffusion == 'MC_2D4pint_ext1D' :
                self.param[2] = resize(self.param[2],(nru,nrv))
                self.nrparam = nru * nrv
                self.rowlength = nrv
                self.pargrid = True
                self.param[3] = resize(self.param[3],(nru,nrv))
                self.param[4] = resize(self.param[4],(nru,nrv))
                self.param[5] = resize(self.param[5],(nru,nrv))
                self.param[0] = uval
                self.param[1] = vval
            elif options.diffusion == 'MC_2D4pint_func':
                grid = eval(options.diffparam[4])
                self.umin = float(grid[0][0])
                self.umax = float(grid[0][1]); nru = grid[0][2]
                self.vmin = float(grid[1][0]);
                self.vmax = float(grid[1][1]); nrv = grid[1][2]
                nru_c = complex(nru,0); nrv_c = complex(nrv,0)
                uval = [float(x) for x in mgrid[self.umin:self.umax:nru_c]]
                vval = [float(x) for x in mgrid[self.vmin:self.vmax:nrv_c]]
                self.nrparam = nru * nrv
                self.rowlength = nrv
                self.pargrid = True
                tD11 = eval(options.diffparam[0])
                tD12 = eval(options.diffparam[1])
                tD21 = eval(options.diffparam[2])
                tD22 = eval(options.diffparam[3])
                for u in uval :
                    for v in vval :
                        self.param[2].append(float(tD11(u,v)))
                        self.param[3].append(float(tD12(u,v)))
                        self.param[4].append(float(tD21(u,v)))
                        self.param[5].append(float(tD22(u,v)))
                self.param[2] = resize(self.param[2],(nru,nrv))
                self.param[3] = resize(self.param[3],(nru,nrv))
                self.param[4] = resize(self.param[4],(nru,nrv))
                self.param[5] = resize(self.param[5],(nru,nrv))
                self.param[0] = uval
                self.param[1] = vval
                        
        except IndexError:
            print 'Too few parameters given for diffusion type', options.diffusion
            sys.exit()
            
        #with the parameters stored, we now make the function
        self.make_diff_func()
            
        #test
        #res=11
        #nr = complex(0,res)
        #stepu = (5.0e-5 - 0. )/ (res-1)
        #stepv = (0.8e-5 - 1.8e-5 )/ (res-1)
        #unew,vnew = mgrid[0.:5.0e-5:nr,0.8e-5:1.8e-5:nr]
        #test = self.diff11(unew[:,0],vnew[0,:])
        #print 'test: ' , test
        #for i in range(res) :
        #    for j in range(res) :
        #        if test[i][j] < 0 :
        #            print " diffusion < 0: [u,v]=[",unew[i,0],",", \
        #                vnew[0,j],"], diff = " , test[i][j] 
        #            print self.D(1,1,unew[i,0]+1e-8, vnew[j,0]+1e-8)
        #            sys.exit()
        
    def _extract_par(self, options):
        #we extract the 1D diffusion points for first component
        uval =[float(x) for x in options.diffparam[0].replace('[',' ').replace(']',' ').split(',')]
        tD11 =[float(x) for x in options.diffparam[1].replace('[',' ').replace(']',' ').split(',')]
        tD12 =[float(x) for x in options.diffparam[2].replace('[',' ').replace(']',' ').split(',')]
        vval =[float(x) for x in options.diffparam[3].replace('[',' ').replace(']',' ').split(',')]
        tD21 =[float(x) for x in options.diffparam[4].replace('[',' ').replace(']',' ').split(',')]
        tD22 =[float(x) for x in options.diffparam[5].replace('[',' ').replace(']',' ').split(',')]
        # we extract vmin, vmax, nrv, umin, umax, nru
        self.umin = uval[0]; self.umax = uval[-1]
        self.vmin = vval[0]; self.vmax = vval[-1]
        nru = len(uval); nrv = len(vval)
        return (uval, tD11, tD12, vval, tD21, tD22, nru, nrv)
        
    def make_diff_func(self) :
        ''' everytime parameters are changed, the diffusion 
            function must be REMADE or it's data adapted !! 
        '''
        if  (self.difftype == 'MC_2D4pint_ext1D' 
                    or self.difftype == 'MC_2D4pint_func'):
                from GridUtils import GridFunc2D
                self.diff11   = GridFunc2D([self.param[0],self.param[1]],N.asarray(self.param[2]),False)
                self.diff12   = GridFunc2D([self.param[0],self.param[1]],N.asarray(self.param[3]),False)
                self.diff21   = GridFunc2D([self.param[0],self.param[1]],N.asarray(self.param[4]),False)
                self.diff22   = GridFunc2D([self.param[0],self.param[1]],N.asarray(self.param[5]),False)
                self.diff11_u = None;self.diff12_u = None;self.diff21_u = None;self.diff22_u = None;
                self.diff11_v = None;self.diff12_v = None;self.diff21_v = None;self.diff22_v = None;
                #grideval means a ugrid and vgrid is given, the grid made by this is evaluated
                #this is used in the plotting function
                self.diff11_grideval   = GridFunc2D([self.param[0],self.param[1]],N.asarray(self.param[2]),True)
                self.diff12_grideval   = GridFunc2D([self.param[0],self.param[1]],N.asarray(self.param[3]),True)
                self.diff21_grideval   = GridFunc2D([self.param[0],self.param[1]],N.asarray(self.param[4]),True)
                self.diff22_grideval   = GridFunc2D([self.param[0],self.param[1]],N.asarray(self.param[5]),True)


        elif self.difftype == 'MC_bspline_ext1D' or \
                  self.difftype == 'MC_bspline1storder_ext1D' : 
                self.arrayaware = False           
                weight = ones(self.nrparam)
                # standard dev on diffusion = 0.2, inverse is 5
                weight[:] =5 * weight[:] 
                #create an interp object
                from scipy.interpolate.fitpack import bisplrep,bisplev
                #param[0] should be list of u values
                #param[1] should be list of v values
                #param[2] should be list of D(u,v) values
                #create the data of the bicubic bspline:
                # no smoother so s = 0
                if self.difftype == 'MC_bspline1storder_ext1D' :
                    self.splinedataD11 = bisplrep(self.param[0], self.param[1],
                        self.param[2], w= weight, s=None,
                        kx = 1, ky=1, full_output = 1, quiet = 0)
                else :
                    #nknots = 11
                    #tx = range(nknots)
                    #tx = [self.umin + x/(1.*(nknots-1)) * (self.umax-self.umin) for x in tx]
                    #ty = range(nknots)
                    #ty = [self.vmin + x/(1.*(nknots-1)) * (self.vmax-self.vmin) for x in ty]
                    #tx = tx[0]*3 + tx + tx[-1]*3
                    #tx[:1] = [tx[0]]*4 ; tx[-1:] = [tx[-1]]*4
                    #ty[:1] = [ty[0]]*4 ; ty[-1:] = [ty[-1]]*4
                    #print 'tx,ty',tx, ty
                    self.splinedataD11 = bisplrep(self.param[0], self.param[1],
                        self.param[2], w= weight, s=2100, 
                        kx = 3, ky=3, full_output = 1, quiet = 0)
                #print 'D11 output :', self.splinedataD11[1:]
                if self.splinedataD11[2] > 0 :
                    print 'ERROR in spline interpolation'
                    sys.exit()
                self.splinedataD11 = self.splinedataD11[0]
                self.diff11   = lambda u,v: bisplev([u], [v], self.splinedataD11, \
                            dx = 0, dy=0)
                self.diff11_u = lambda u,v: bisplev([u], [v] , self.splinedataD11,\
                            dx = 1, dy=0)
                self.diff11_v = lambda u,v: bisplev([u], [v] , self.splinedataD11,\
                            dx = 0, dy=1)
                self.diff11_grideval   = lambda u,v: bisplev(u, v, self.splinedataD11, \
                            dx = 0, dy=0)
                self.diff11_u_grideval = lambda u,v: bisplev(u, v , self.splinedataD11,\
                            dx = 1, dy=0)
                self.diff11_v_grideval = lambda u,v: bisplev(u, v , self.splinedataD11,\
                            dx = 0, dy=1)
                if self.difftype == 'MC_bspline1storder_ext1D' :
                    self.splinedataD12 = bisplrep(self.param[0], self.param[1],
                        self.param[3], w= weight, s=None,
                        kx = 1, ky=1, full_output = 1, quiet = 1)
                else :
                    self.splinedataD12 = bisplrep(self.param[0], self.param[1],
                        self.param[3], w= weight, s=None,
                        kx = 3, ky=3, full_output = 1, quiet = 0)
                #print 'D12 output :', self.splinedataD12[1:]
                if self.splinedataD12[2] > 0 :
                    print 'ERROR in spline interpolation'
                    sys.exit()
                self.splinedataD12 = self.splinedataD12[0]
                self.diff12   = lambda u,v: bisplev([u], [v], self.splinedataD12, \
                            dx = 0, dy=0)
                self.diff12_u = lambda u,v: bisplev([u], [v] , self.splinedataD12,\
                            dx = 1, dy=0)
                self.diff12_v = lambda u,v: bisplev([u], [v] , self.splinedataD12,\
                            dx = 0, dy=1)
                self.diff12_grideval   = lambda u,v: bisplev(u, v, self.splinedataD12, \
                            dx = 0, dy=0)
                self.diff12_u_grideval = lambda u,v: bisplev(u, v , self.splinedataD12,\
                            dx = 1, dy=0)
                self.diff12_v_grideval = lambda u,v: bisplev(u, v , self.splinedataD12,\
                            dx = 0, dy=1)
                if self.difftype == 'MC_bspline1storder_ext1D' :
                    self.splinedataD21 = bisplrep(self.param[0], self.param[1],
                        self.param[4], w= weight, s=None,
                        kx = 1, ky=1, full_output = 1, quiet = 1)
                else :
                    self.splinedataD21 = bisplrep(self.param[0], self.param[1],
                        self.param[4], w= weight, s=None,
                        kx = 3, ky=3, full_output = 1, quiet = 0)
                #print 'D21 output :', self.splinedataD21[1:]
                if self.splinedataD21[2] > 0 :
                    print 'ERROR in spline interpolation'
                    sys.exit()
                self.splinedataD21 = self.splinedataD21[0]
                self.diff21   = lambda u,v: bisplev([u], [v], self.splinedataD21, \
                            dx = 0, dy=0)
                self.diff21_u = lambda u,v: bisplev([u], [v] , self.splinedataD21,\
                            dx = 1, dy=0)
                self.diff21_v = lambda u,v: bisplev([u], [v] , self.splinedataD21,\
                            dx = 0, dy=1)
                self.diff21_grideval   = lambda u,v: bisplev(u, v, self.splinedataD21, \
                            dx = 0, dy=0)
                self.diff21_u_grideval = lambda u,v: bisplev(u, v , self.splinedataD21,\
                            dx = 1, dy=0)
                self.diff21_v_grideval = lambda u,v: bisplev(u, v , self.splinedataD21,\
                            dx = 0, dy=1)
                if self.difftype == 'MC_bspline1storder_ext1D' :
                    self.splinedataD22 = bisplrep(self.param[0], self.param[1],
                        self.param[5], w= weight, s=None,
                        kx = 1, ky=1, full_output = 1, quiet = 1)
                else :
                    self.splinedataD22 = bisplrep(self.param[0], self.param[1],
                        self.param[5], w= weight, s=55000.,
                        kx = 3, ky=3, full_output = 1, quiet = 0)
                #print 'D22 output :', self.splinedataD22[1:]
                if self.splinedataD22[2] > 0 :
                    print 'ERROR in spline interpolation'
                    sys.exit()
                self.splinedataD22 = self.splinedataD22[0]
                self.diff22   = lambda u,v: bisplev([u], [v], self.splinedataD22, \
                            dx = 0, dy=0)
                self.diff22_u = lambda u,v: bisplev([u], [v] , self.splinedataD22,\
                            dx = 1, dy=0)
                self.diff22_v = lambda u,v: bisplev([u], [v] , self.splinedataD22,\
                            dx = 0, dy=1)
                self.diff22_grideval   = lambda u,v: bisplev(u, v, self.splinedataD22, \
                            dx = 0, dy=0)
                self.diff22_u_grideval = lambda u,v: bisplev(u, v , self.splinedataD22,\
                            dx = 1, dy=0)
                self.diff22_v_grideval = lambda u,v: bisplev(u, v , self.splinedataD22,\
                            dx = 0, dy=1)
        elif self.difftype == 'MC_interp2d_ext1D':
                from scipy.interpolate.interpolate import interp2d
                self.diff11   = interp2d(self.param[0],self.param[1],self.param[2])
                self.diff12   = interp2d(self.param[0],self.param[1],self.param[3])
                self.diff21   = interp2d(self.param[0],self.param[1],self.param[4])
                self.diff22   = interp2d(self.param[0],self.param[1],self.param[5])
                print 'diff11 test ' , self.diff11(1.5)
                print ' *** DO NOT USE THIS INTERPOLATION, it performs badly ***'
                sys.exit()
        elif self.difftype == 'function' :
            #self.diffij is defined when reading the init file for scalars,
            pass
        else:
            print ('Unknown diffusion type given', self.difftype)
            sys.exit()
        
        #store the functions in a list:
        self.diff = [[self.diff11 , self.diff12],[self.diff21 , self.diff22]]
        self.diff_u = [[self.diff11_u , self.diff12_u],[self.diff21_u , self.diff22_u]]
        self.diff_v = [[self.diff11_v , self.diff12_v],[self.diff21_v , self.diff22_v]]
        #changed parameters are taken into account, so diffusion is now correct
        self.modified = False
            
    def nrParam_uu(self) :
        return self.nrparam
    def nrParam_uv(self) :
        return self.nrparam
    def nrParam_vu(self) :
        return self.nrparam
    def nrParam_vv(self) :
        return self.nrparam
    def get_param_uu(self,number):
        return self._get_param(2,number)
    def get_param_uv(self,number):
        return self._get_param(3,number)
    def get_param_vu(self,number):
        return self._get_param(4,number)
    def get_param_vv(self,number):
        return self._get_param(5,number)
    def _get_param(self, order, number) :
        if self.pargrid :
            return self.param[order][number/self.rowlength][number - (number/self.rowlength)*self.rowlength]
        else :
            return self.param[order][number]
    def get_param_pos(self, number) :
        '''For all orders (uu,uv,vu,vv) the parameters are on the same u-v grid
           This returns a u_pos, vpos, that is the position of the parameter in uv space
        '''
        if self.pargrid :
            return self.param[0][number/self.rowlength],self.param[1][number - (number/self.rowlength)*self.rowlength]
        else :
            return self.param[0][number],self.param[1][number]
    def set_param_uu(self,number,val):
        self._set_param(2, number, val)
    def set_param_uv(self,number,val):
        self._set_param(3, number, val)
    def set_param_vu(self,number,val):
        self._set_param(4, number, val)
    def set_param_vv(self,number,val):
        self._set_param(5, number, val)
    def _set_param(self, order, number, val) :
        if self.pargrid :
            self.param[order][number/self.rowlength][number - (number/self.rowlength)*self.rowlength] = val
            self.modified = True
        else :
            self.param[order][number] = val
            self.modified = True
        
                
    def D(self, i, j, u, v) :
        ''' call interdiffusion coefficient diffij in numbers u, v
        '''
        if self.modified : 
            self.make_diff_func()
        return self.diff[i-1][j-1](u,v)
    
    def D_a(self, i, j, u, v) :
        ''' call interdiffusion coefficient diffij in arrays u, v as u[i],v[i]
        '''
        if self.modified : 
            self.make_diff_func()
        if self.arrayaware :
            return self.diff[i-1][j-1](u,v)
        else :
            sol = [0.]*len(u)
            for k in range(len(sol)):
                sol[k]=self.diff[i-1][j-1](u[k],v[k])
            return sol
    
    def D_u(self, i, j, u, v) :
        ''' call deriv interdiffusion coefficient diffij_u in numbers u, v
        '''
        if self.modified : 
            self.make_diff_func()
        return self.diff_u[i-1][j-1](u,v)
    
    def D_u_a(self, i, j, u, v) :
        ''' call deriv interdiffusion coefficient diffij_u in arrays u, v
        '''
        if self.modified : 
            self.make_diff_func()
        if self.arrayaware :
            return self.diff_u[i-1][j-1](u,v)
        else :
            sol = [0]*len(u)
            for k in range(len(sol)):
                sol[k]=self.diff_u[i-1][j-1](u[k],v[k])
            return sol
    
    def D_v(self, i, j, u, v) :
        ''' call deriv interdiffusion coefficient diffij_v in numbers u, v
        '''
        if self.modified : 
            self.make_diff_func()
        return self.diff_v[i-1][j-1](u,v)
    
    def D_v_a(self, i, j, u, v) :
        ''' call deriv interdiffusion coefficient diffij_v in arrays u, v
        '''
        if self.modified : 
            self.make_diff_func()
        if self.arrayaware :
            return self.diff_v[i-1][j-1](u,v)
        else :
            sol = [0]*len(u)
            for k in range(len(sol)):
                sol[k]=self.diff_v[i-1][j-1](u[k],v[k])
            return sol        
    
    def plot(self, res=70, show=False, dir='test', ext='0'):
        ''' plot D(u,v) and diff(D,u) between min and max boundaries
            and save the plot to a .vtk file. 
            It saves to path.join(dir,dir+'MCDiff'+ext+'.vtk') 
            Open file with mayavi or paraview
            Resolution of 70
        '''
        if self.modified : 
            self.make_diff_func()
        if not os.path.isdir(dir):
            print ("Error, plot diffusion not possible, directory does not exist")
            sys.exit()
        
        #nrpoints
        nr = complex(0,res)
        
        #we make equal x,y,z axes for vtk plot:
        stepu = (self.umax - self.umin )/ (res-1)
        stepv = (self.vmax - self.vmin )/ (res-1)
        unew,vnew = mgrid[self.umin:self.umax:nr,self.vmin:self.vmax:nr]
        
        if self.difftype == 'function' :
            #create the arrays to hold gridfunction
            from GridUtils import GridFunc2D
            uval = list(mgrid[self.umin:self.umax:nr])
            vval = list(mgrid[self.vmin:self.vmax:nr])
            listu = []
            listv = []
            for i in range(0,res) :
                    for j in range(0,res) :
                        listu.append(uval[i])
                        listv.append(vval[j])
            val = self.D_a(1,1, listu, listv)
            self.diff11_grideval = GridFunc2D([uval,vval],N.asarray(val),True)
            val = self.D_a(1,2, listu, listv)
            self.diff12_grideval = GridFunc2D([uval,vval],N.asarray(val),True)
            val = self.D_a(2,1, listu, listv)
            self.diff21_grideval = GridFunc2D([uval,vval],N.asarray(val),True)
            val = self.D_a(2,2, listu, listv)
            self.diff22_grideval = GridFunc2D([uval,vval],N.asarray(val),True)
            
        #plot with gnuplot
        g  = Gnuplot.Gnuplot(persist=0)
        g.title('D11 Interdiffusion coefficient')
        g.xlabel('Al conc [mol/mm^3]')
        g.ylabel('Si conc  [mol/mm^3]')
        g.__call__('set zlabel "D11')
        g.__call__('set xtics 1e-5')
        g.__call__('set ytics 1e-5')
        g.__call__('set ztics 2')
        #print self.diff11_grideval(unew[:,0],vnew[0,:])
        g.splot(Gnuplot.GridData(self.diff11_grideval(unew[:,0],vnew[0,:]),
                    unew[:,0],vnew[0,:],binary=True))
        g.hardcopy(filename=os.path.join(dir,dir+'_D11.ps'), enhanced=1,mode='eps',
                color=0,fontname='Times-Roman',fontsize=28)
        g.title('D12 Interdiffusion coefficient')
        g.__call__('set zlabel "D12')
        g.splot(Gnuplot.GridData(self.diff12_grideval(unew[:,0],vnew[0,:]),
                    unew[:,0],vnew[0,:],binary=True))
        g.hardcopy(filename=os.path.join(dir,dir+'_D12.ps'), enhanced=1,mode='eps',
                color=0,fontname='Times-Roman',fontsize=28)
        g.title('D21 Interdiffusion coefficient')
        g.__call__('set zlabel "D21')
        g.splot(Gnuplot.GridData(self.diff21_grideval(unew[:,0],vnew[0,:]),
                    unew[:,0],vnew[0,:],binary=True))
        g.hardcopy(filename=os.path.join(dir,dir+'_D21.ps'), enhanced=1,mode='eps',
                color=0,fontname='Times-Roman',fontsize=28)
        g.title('D22 Interdiffusion coefficient')
        g.__call__('set zlabel "D22')
        g.splot(Gnuplot.GridData(self.diff22_grideval(unew[:,0],vnew[0,:]),
                    unew[:,0],vnew[0,:],binary=True))
        g.hardcopy(filename=os.path.join(dir,dir+'_D22.ps'), enhanced=1,mode='eps',
                color=0,fontname='Times-Roman',fontsize=28)
    
        znew11 = self.diff11_grideval(unew[:,0],vnew[0,:])
        znew12 = self.diff12_grideval(unew[:,0],vnew[0,:])
        znew21 = self.diff21_grideval(unew[:,0],vnew[0,:])
        znew22 = self.diff22_grideval(unew[:,0],vnew[0,:])
        # Flatten the 2D array data as per VTK's requirements.
        z11 = reshape(transpose(znew11), (-1,))
        z12 = reshape(transpose(znew12), (-1,))
        z21 = reshape(transpose(znew21), (-1,))
        z22 = reshape(transpose(znew22), (-1,))
        # now dump the data to a VTK file.
        import pyvtk
        s11 = pyvtk.Scalars(z11, name = '11', lookup_table = 'default')
        s21 = pyvtk.Scalars(z12, name = '21', lookup_table = 'default')
        s12 = pyvtk.Scalars(z21, name = '12', lookup_table = 'default')
        s22 = pyvtk.Scalars(z22, name = '22', lookup_table = 'default')
        s11.default_float = 'double'; s12.default_float = 'double'
        s21.default_float = 'double'; s22.default_float = 'double'
        point_data11 = pyvtk.PointData(s11)
        point_data12 = pyvtk.PointData(s12)
        point_data21 = pyvtk.PointData(s21)
        point_data22 = pyvtk.PointData(s22)
        
        #StructuredPoints: dimension of grid (nr points in x, y, z dir), 
        #                  origin of made grid (we take it in 0)
        #                  spacing of grid
        # We do not need z axis so we obtain a 2d grid, with z component = 0
        # We scale the grid to have a better picture
        scale = 1e6
        grid = pyvtk.StructuredPoints((res,res, 1), #(self.umin*scale, self.vmin*scale, 0),\
                                      (0.,0.,0.),
                                      (stepu*scale, stepv*scale, 10e-6))
        grid.default_float = 'double'
        #Data over this grid, these are scalars given in the points
        data11 = pyvtk.VtkData(grid, 'Interdiffusion Al-Si D11', point_data11)
        data12 = pyvtk.VtkData(grid, 'Interdiffusion Al-Si D12', point_data12)
        data21 = pyvtk.VtkData(grid, 'Interdiffusion Al-Si D21', point_data21)
        data22 = pyvtk.VtkData(grid, 'Interdiffusion Al-Si D22', point_data22)
        
        data11.tofile(os.path.join(dir,dir+'MCDiff11'+ext+'.vtk'))
        data12.tofile(os.path.join(dir,dir+'MCDiff12'+ext+'.vtk'))
        data21.tofile(os.path.join(dir,dir+'MCDiff21'+ext+'.vtk'))
        data22.tofile(os.path.join(dir,dir+'MCDiff22'+ext+'.vtk'))
        
        '''
        showMayavi = False
        if showMayavi :
            import mayavi
            v = mayavi.mayavi() # create a MayaVi window.
            if show :
                d = v.open_vtk(os.path.join(dir,dir+'MCDiff'+ext+'.vtk', config=0)) # open the data file.
            else :
                d = v.open_vtk(os.path.join(dir,dir+'MCDiff'+ext+'.vtk', config=1)) # open the data file.
            # The config option turns on/off showing a GUI control for the data/filter/module.
            # load the filters.
            f = v.load_filter('WarpScalar', config=0) 
            n = v.load_filter('PolyDataNormals', 0)
            n.fil.SetFeatureAngle (45) # configure the normals.
            # Load the necessary modules.
            m = v.load_module('SurfaceMap', 0)
            a = v.load_module('Axes', 0)
            a.axes.SetCornerOffset(0.0) # configure the axes module.
            o = v.load_module('Outline', 0)
            v.Render() # Re-render the scene.
        else :
            print 'file', os.path.join(dir,dir+'MCDiff'+'.vtk'), \
                    'written. Open it to see diffusion' 
        '''
    

class Diffusion_xt_uxt_vxt:
    ''' a diffusion function, returning given the time t and position x, 
        the diffusion D(x,t), using a function u and v, so internally D(u,v) 
        can be calculated, knowing u(x,t) and v(x,t) gives D(u(x,t), v(x,t)) 
    '''
    def __init__(self, gridx, deltax, griduval = True, equidistant=True):
        '''
        set the diffusion over a grid x, which is a numpy list.
        The list should be equidistant, then deltax is used
        If griduval = True : data of u,v is in the gridpoints
                     = False: data of u,v is in the cell midpoint
        '''
        self.grid = gridx
        self.lengrid = len(gridx)
        self.stepgrid = deltax
        self.diff = None
        self.griduval = griduval
        self.uval_time = None
        #only equidistant implemented
        assert(equidistant)
        #only cell data implemented
        assert(not(griduval))
        
        
    def set_diff_func(self, diff):
        ''' Set the diffusion function in u, v used, eg a Diffusion_u_MC class 
            object, should have components 11, 12, 21, 22
        '''
        self.diff = diff
        
    def set_grid_func_time_method(self, methodu_v):
        ''' set the method with wich uxt, vxt at a required time can be obtained,
            so method(time) returns an uval
        '''
        self.uv_val_time_method = methodu_v
        
    def set_grid_func(self, uval, vval, time, griduval=True):
        ''' diffusion will be returned based on the values uval, vval; 
            which are values correct for time
            If griduval, the given uval/vval is in gridpoints, 
            if False it is in cell
        '''
        if self.griduval == griduval:
            assert(len(uval) == len(self.grid))
            assert(len(vval) == len(self.grid))
        else :
            assert(len(uval) == len(self.grid)-1)
            assert(len(vval) == len(self.grid)-1)
        self.uval = uval
        self.vval = vval
        self.uval_time = time
        if self.griduval == griduval:
            #the uval values corresponding to self.grid
            self.uval_grid = uval
            self.vval_grid = vval
        else :
            print 'diffusion.py-set_grid_func, not implemented '
            sys.exit()
        
    def D(self, time, xval, i, j):
        ''' return value of diffusion is (numpy) array of points xval
            for diff component ij
        '''
        #determine the x coord needed, following are normal slow arrays as result
        if not(time == self.uval_time) :
            #we need to get uval at time 
            uval, vval = self.uv_val_time_method(time)
            self.set_grid_func(uval, vval, time, self.griduval)
            
        def indexcorrect(index) :
            i=int(index)
            if i == self.lengrid-1 :
                i -= 1
            return i
        
        #we use the fact the grid is equidistant
        xpos = map(indexcorrect,(xval - self.grid[0])/self.stepgrid)
        uval = array(map(lambda x: self.uval_grid[x], xpos))
        uval_right = array(map(lambda x: self.uval_grid[x+1], xpos))
        vval = array(map(lambda x: self.vval_grid[x], xpos))
        vval_right = array(map(lambda x: self.vval_grid[x+1], xpos))
        grid = array(map(lambda x: self.grid[x], xpos))
        grid_right = array(map(lambda x: self.grid[x+1], xpos))
        #return interpolated value of D in x, so u is with (x_i,u_i) and  (x_i+1,u_i+1)
        uvalgrid = (uval_right -uval)/(grid_right-grid)*(xval-grid) + uval
        vvalgrid = (vval_right -vval)/(grid_right-grid)*(xval-grid) + vval
            
        return self.diff.D_a(i, j, uvalgrid, vvalgrid)
