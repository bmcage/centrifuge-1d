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

import sys
import numpy as N
from scipy import *
from scipy import weave
from scipy.weave import converters
from bisect import bisect_left
import copy

#-------------------------------------------------------------------------
#
# Some Constants for grids
#
#-------------------------------------------------------------------------

def DL_1d(grid,full=True):
    ''' 1d grid first derivative over a grid
            grid : a list
    '''
    size = len(grid)
    if full :
        mat = sparse.lil_matrix((size,size)) 
        mat[0,0]   = -1/(grid[1]-grid[0])
        mat[0,1]   =  1/(grid[1]-grid[0])
        for i in range(1,size-1,1) :
            mat[i,i-1] = -1/(grid[i+1]-grid[i-1])
            mat[i,i+1] =  1/(grid[i+1]-grid[i-1])
        mat[size-1,size-2]   = -1/(grid[-1]-grid[-2])
        mat[size-1,size-1]   =  1/(grid[-1]-grid[-2])
    else :
        mat = sparse.lil_matrix((size-1,size)) 
        for i in range(size-1) :
            mat[i,i]   = -1/(grid[i+1]-grid[i])
            mat[i,i+1] =  1/(grid[i+1]-grid[i])
    return mat.tocsc()


#-------------------------------------------------------------------------
#
# Some Grid functions
#
#-------------------------------------------------------------------------

class DomainError(Exception):
    def __init__(self,x):
        self.x = x
    def __str__(self):
        if len(self.x) == 3 :
            return 'x=%g is not in [%g,%g]' %(x[0],x[1],x[2])
        else :
            return 'Data given not inside of domain'
#-------------------------------------------------------------------------
#
# Create a 1D grid, 
#
#-------------------------------------------------------------------------

def create_grid(gridpoints=11, begin=0., end=1., type="equidistant",\
                refine=1.):
    """Make a grid from begin to end consisting of gridpoints points
       if type = equidistant : equal interval length
       if type = borderrefined: at borders refine times smaller than center    
    """
    step = (end-begin)/(gridpoints-1)
    if type == "equidistant" or (type=="borderrefined" and refine==1.):
        g, step = linspace(begin, end, num=gridpoints, retstep=True)
        return g;
    elif type == "borderrefined":
        #power will be such that in the symmetric case we almost have that
        #power * shortest interval = longest interval dense points at begin and end
        pw = 2 * math.log(refine) /( gridpoints-1);
        g = zeros(gridpoints, float)
        for i in xrange(g.shape[0]):
            if i==0: continue
            j = gridpoints-i
            if j> i : j=i
            #we do previous value plus interval length. This length is such:
            #for i=1=N we have length1= power^(2/N) = 1 (large N), 
            #for i=N/2: length2=power, so longest interval is power*shortest
            g[i]=g[i-1]+math.exp(pw * j)
        lastval = g[gridpoints-1]
        #normalize distribution on required interval
        g = begin + (end - begin) * g / lastval;
        return g
    else:
        print "invalid option in create grid, type=",type
        sys.exit()

def pos_in_grid(posrel_0_1, cutleft, length):
    """posrel is the relative position in [0,1], length is the real length 
       of the interval, and cutleft the cutoff to the left"""
    return cutleft + posrel_0_1*length

#-------------------------------------------------------------------------
#
# Project a concentration to a fixed grid
#
#-------------------------------------------------------------------------

def read_exp_file(filename):
    """
    Read in experiment concentration file
    """
    infile = open(filename,'r')
    line = infile.readline().split()
    time = float(line[0])
    nrlines  = int(line[1])
    temp     = float(line[2])
    wtp_Si   = float(line[3])
    line = infile.readline()  # dummy line
    meshinitpos = []               
    concAlinit  = []
    concSiinit  = []
    while 1:                  # read in the file
        line = infile.readline()
        if not line: break
        pos, concAl, concSi =line.split()
        meshinitpos.append(float(pos))  #in micrometer
        concAlinit.append(float(concAl))
        concSiinit.append(float(concSi))
        
    infile.close()
    if len(meshinitpos) != nrlines:
        print ("inputfile '"+ filename +"' contains wrong number of lines")
    return time, temp, wtp_Si, meshinitpos, concAlinit, concSiinit
    
#-------------------------------------------------------------------------
#
# Project a concentration to a fixed grid
#
#-------------------------------------------------------------------------

def project_exp_to_grid(meshexp, valexp, grid,grid_cutoff,grid_length):
    '''
    project a mesh meshexp with experiment values valexp onto grid
    Here grid is a grid between [0,1], corresponding to 
        [grid_cutoff, grid_cutoff+grid_length]
    Return value: a numarray with the projected values, contains None for
                  grid values outside of experiment data range
    '''
    # the storage
    nrp = len(grid)
    valexpongrid = zeros(nrp, float)
    # the exp sample might be smaller than the grid, set value None there
    b=0
    while meshexp[0]>grid[b]:
        valexpongrid[b] = None
        b += 1
    e=nrp-1
    while meshexp[-1]<grid[e]:
        valexpongrid[e] = None
        e -= 1
    j=1
    for i in range(b,e,1):
        valexpongrid[i] = valexp[j-1] + (valexp[j]-valexp[j-1]) \
            /(meshexp[j]-meshexp[j-1]) \
            * (pos_in_grid(grid[i],grid_cutoff,grid_length)-meshexp[j-1])
        while pos_in_grid(grid[i+1],grid_cutoff,grid_length) > meshexp[j]:
            j+=1
            if j == len(meshexp):
                break
    #last val
    if i != e-1 :
        print ("ERROR, projection measurement data on grid failed.")
        print ("i is",i, "e is ",e)
        sys.exit()
    i = e
    if j == len(meshexp):
        #check if ok and set j correct
        j = j-1
        if not (  meshexp[j-1] <= pos_in_grid(grid[i],grid_cutoff,grid_length) and 
                  (1+eps)*meshexp[j] > pos_in_grid(grid[i],grid_cutoff,grid_length) ):
            print ("ERROR, projection measurement data on grid failed.")
            print ("check", meshexp[j-1],'<=',pos_in_grid(grid[i],grid_cutoff,grid_length),'<=',meshexp[j])
            sys.exit()
    valexpongrid[i]= valexp[j-1] + (valexp[j]-valexp[j-1]) \
                        /(meshexp[j]-meshexp[j-1]) \
                        * (pos_in_grid(grid[i],grid_cutoff,grid_length)-meshexp[j-1])
                        
    return valexpongrid

#-------------------------------------------------------------------------
#
# fuse 2 solutions over a grid to one array
#
#-------------------------------------------------------------------------
def fuse(u,v) :
    if len(u) != len(v) :
        print 'ERROR in GridUtils.fuse '
        sys.exit()
    uv = zeros(2*len(u), float)
    for i in range(len(u)) :
        uv[2*i]   = u[i]
        uv[2*i+1] = v[i]
    return uv

#-------------------------------------------------------------------------
#
# defuse 2 fused solutions over a grid to a tupple of 2 solutions
#
#-------------------------------------------------------------------------
def defuse(uv) :
    if mod(len(uv),2) != 0 :
        print 'ERROR in GridUtils.defuse'
        sys.exit()
    u = zeros(len(uv)/2,float)
    v = zeros(len(uv)/2,float)
    for i in range(len(uv)/2) :
        u[i] = uv[2*i]
        v[i] = uv[2*i+1]
    return u , v

#-------------------------------------------------------------------------
#
# A class for functions over a grid
#
#-------------------------------------------------------------------------

class GridFunc:
    ''' A class to hold interpolation function of a value field over a rectangular grid.
        Can be called with data in the grid
    '''
    def __init__(self,dim=1, grid=None, val=None):
        ''' init of a gridfunction. The dim must be given, and grid and val
            must be arrays. grid an array[0..dim-1] with the grid values on coordinate axis
            and val an array with a matrix of the values, each index the index of the 
            corresponding coordinate axis point.
        '''
        self.dim=dim
        self.grid=grid
        self.val=val
        self.gridsize_x = len(self.grid[0])
        
    def grid(self):
        return self.grid
    def values(self):
        return self.val
    def dim(self):
        return self.dim

class GridFunc1D(GridFunc):
    ''' A class to hold interpolation function of a value field over a 1D grid.
        Can be called with data in the grid
    '''
    def __init__(self,grid=None, val=None):
        '''
        *grid is a list, containing one list
        *val is a list containing function values
        '''
        if len(grid)==1 and len(val)==len(grid[0]) :
            GridFunc.__init__(self,1,grid,val)
        else : 
            print 'ERROR: 1D gridfunc with wrong data,',len(grid), len(val), len(grid[0])
            raise ValueError
            sys.exit()
        self.lastxpos = 0
            
    def __call__(self,pos=None):
        '''call of the function with in position pos, returns the value
        '''
        if isinstance(pos,(N.ndarray, ndarray, list)) : 
            #return [self._call_val_weave(x) for x in pos]
            return self._call_val_weave(N.asarray(pos)) 
        else :
            return self._call_val_bisect(pos)
        
    def _domain_pos(self,pos):
        '''Return the index so that pos is between index and index+1
           Raise DomainError if not in the domain
        '''
        if pos<self.grid[0][0] or pos > self.grid[0][-1] :
            raise DomainError, [pos, self.grid[0][0], self.grid[0][-1]]
        
        xpos = self.lastxpos
        if self.grid[0][xpos] <= pos and self.grid[0][xpos+1]>= pos :
            return xpos
        xpos = bisect_left(self.grid[0],pos)
        xpos -= 1
        if xpos == -1 : xpos = 0
        return xpos
        
    def map(self, pos):
        ''' Given a pos between begin and end, give mapping to values in 
             nodes:
              gridfunc(x) = sum(gridval_i * phi_i(x)) where phi_i is localised wavelet in point grid_i
            Generally phi_i(x) = 0 exept for two points, say x_k and x_k+1
            This function returns the map (x_k,phi_i(x)) as a sparse tuple, that is
              nr_of_returnval, listindex, mapval, eg:
                (2 , [4, 5] , [0.7, 0.3] ) 
                  --> this indicates that in grid[4] we have 0.7 * value[4] and in grid[5] we have 0.3*value[5]
        '''
        xpos = self._domain_pos(pos)
        self.lastxpos = xpos
        #return the value
        return (2 , [xpos, xpos+1], 
             [(self.grid[0][xpos+1]-pos)/(self.grid[0][xpos+1]-self.grid[0][xpos])
               , (pos-self.grid[0][xpos])/(self.grid[0][xpos+1]-self.grid[0][xpos]) ]
             )
      
    def _call_val_bisect(self,pos):
        #first find the index where the data is
        if pos<self.grid[0][0] or pos > self.grid[0][-1] :
            print 'ERROR: value of func asked in ', pos, 'outside of domain [', self.grid[0][0],',',self.grid[0][-1],']'
            sys.exit()
        xpos = self.lastxpos
        if self.grid[0][xpos] <= pos and self.grid[0][xpos+1]>= pos :
            return self.val[xpos] + (self.val[xpos+1]-self.val[xpos])/(self.grid[0][xpos+1]-self.grid[0][xpos]) \
                                    *(pos-self.grid[0][xpos])
        xpos = bisect_left(self.grid[0],pos)
        xpos -= 1
        if xpos == -1 : xpos = 0
        self.lastxpos = xpos
        #interpolate
        return self.val[xpos] + (self.val[xpos+1]-self.val[xpos])/(self.grid[0][xpos+1]-self.grid[0][xpos]) \
                                    *(pos-self.grid[0][xpos])
                                      
    def _call_val(self,pos):
        #first find the index where the data is
        if pos<self.grid[0][0] or pos > self.grid[0][-1] :
            print 'ERROR: value of func asked in ', pos, 'outside of domain [', self.grid[0][0],',',self.grid[0][-1],']'
            sys.exit()
        xpos = self.lastxpos
        if self.grid[0][xpos] <= pos and self.grid[0][xpos+1]>= pos :
            return self.val[xpos] + (self.val[xpos+1]-self.val[xpos])/(self.grid[0][xpos+1]-self.grid[0][xpos]) \
                                    *(pos-self.grid[0][xpos])
        xpos=0
        for gridx in self.grid[0][1:] :
            if pos<=gridx:
                break
            xpos +=1
        self.lastxpos = xpos
        #interpolate
        return self.val[xpos] + (self.val[xpos+1]-self.val[xpos])/(self.grid[0][xpos+1]-self.grid[0][xpos]) \
                                    *(pos-self.grid[0][xpos])
                                    
    def _call_val_weave(self,pos):
        # input is an numpy array
        #   No check on data integrity here
        xvalarray = pos
        nr   = len(pos)
        result = zeros(nr, float)
        xvals= self.grid[0]
        nxval= len(self.grid[0])
        xvalmin = self.grid[0][0]
        xvalmax = self.grid[0][-1] 
        yvals= self.val
        assert(type(xvals)==type([]))
        code = r"""
#line 297 "_call_val_weave"
double xvalsi,yvalsi,xvalsip, yvalsip,xval;
int index,tmp;
int ok;
index=0;
tmp=0;
ok=1;
for (int j=0; j<nr; j++) {
    xval = xvalarray[j];
    if (xval < xvalmin || xval > xvalmax ) {
        ok=0;
        j=nr;
        for (int k=0; k<nr;k++) {
            result[k]=0.0001;
            }
        } else 
    {
    if (xval>= py_to_float(PyList_GetItem(xvals,tmp),"xvalt1") && xval<= py_to_float(PyList_GetItem(xvals,tmp+1),"xvalt2")) {
        index=tmp;
        } else if (tmp+1 < nxval-1 &&  xval>= py_to_float(PyList_GetItem(xvals,tmp+1),"xvalt3") && xval<= py_to_float(PyList_GetItem(xvals,tmp+2),"xvalt4")) {
        index=tmp+1;
        tmp=index;
        } else {
            for (int i=1; i<nxval; i++){
              if (xval <= py_to_float(PyList_GetItem(xvals,i),"xvalsip")) {
                  index = i-1;
                  tmp = index;
                  i = nxval;
                  }
              }
        }
    xvalsip= py_to_float(PyList_GetItem(xvals,index+1),"xvalsip");
    xvalsi = py_to_float(PyList_GetItem(xvals,index),"xvalsi");
    yvalsi = py_to_float(PyList_GetItem(yvals,index),"yvalsi");
    yvalsip = py_to_float(PyList_GetItem(yvals,index+1),"yvalsip");
    result[j] = yvalsi + (yvalsip-yvalsi)/(xvalsip-xvalsi)*(xval-xvalsi);
    }
}
return_val = ok;
"""
        # compiler keyword only needed on windows with MSVC installed
        msg = weave.inline_tools.inline(code,
                           ['nr', 'xvalarray', 'nxval', 'xvals', 'yvals','result','xvalmin','xvalmax']
                           #,type_converters = converters.blitz
                           ,verbose=0)
        if not msg == 1 :
            print 'ERROR error in _call_val_weave, out of bounds' +\
                  ' bounds: ',xvalmin,xvalmax,'data:', xvalarray, \
                  ' In case of optimization: reduce stepsize.'
            raise DomainError, []
            #sys.exit(0)
        return result
                                        
class GridFunc2D(GridFunc):
    ''' A class to hold interpolation function of a value field over a 2D rectangular grid.
        Can be called with data in the grid.
    '''
    def __init__(self,grid, val,gridcall=True):
        '''grid is [xdiv, ydiv], with xdiv array containing points on x axis, idem ydiv
           val is array of size = size(xdiv) x size(ydiv) containing all values in matrix
           *grid is a list of lists!
           *val is a numpy array with two dimensions
           *gridcall: if True, calling the func is taken to be grid evaluation: [1,2],[1,2] evals 11,12,21,22,
                                 which will be returned as a matrix!
                      if False, calling the func is taken to be a list evaluation: [1,2],[1,2] evals 11,22,
                                 which will be returned as a list.
        '''
        if len(grid)==2  and size(val)==size(grid[0])*size(grid[1]):
            GridFunc.__init__(self,2,grid,val)
        else : 
            print 'ERROR GRIDUTILS, len grid', len(grid), 'val',size(val),size(grid[0]),size(grid[1])
            print 'ERROR: 2D gridfunc with wrong data'
            
        self.lenx = len(grid[0])
        self.leny = len(grid[1])
        self.lastxpos = 0
        self.lastypos = 0
        self.gridcall = gridcall
        
        
    def __call__(self,posx, posy, gridcall=None):
        '''call of the function in position (posx,posy), returns the value
           A 4-point interpolation is used:
             phi_1 = 1/4*(1+xa)(1+ya)
             phi_2 = 1/4*(1-xa)(1+ya)
             phi_3 = 1/4*(1-xa)(1-ya)
             phi_4 = 1/4*(1+xa)(1-ya)
           with 
             xa=2 (x-x0)/a
             ya=2 (y-y0)/b
             x0,y0 is central point
             a,b is length and height.
           where the rectangle is numbered
                2 ----- 1
                |       |
                |       |
                3 ----- 4
                
          If posx/y is an array, and grid=True, then the call is done in posx x posy,
          if grid=False, the call is done in all (posx[i],posy[i])
        '''
        if not gridcall == None :
            grid = gridcall
        else :
            grid = self.gridcall
            
        if isinstance(posx,(N.ndarray, ndarray, list)) :
            if grid:
                #print 'called as weave_grid'
                return self._call_val_weave_grid(posx,posy) 
            else :
                #print 'called as weave'
                if not size(posx)==size(posy) :
                    print 'ERROR: input GridFunc2D wrong'
                    sys.exit()
                return self._call_val_weave(posx,posy) 
            '''
            if grid :
                r=zeros(size(posx)*size(posy),'f')
                for i in range(size(posx)) :
                    for j in range(size(posy)) :
                        r[i*size(posx)+j] = self._call_val_bisect(posx[i],posy[j])
                return r
            else :
                n=size(posx)
                if not n==size(posy) :
                    print 'ERROR: input GridFunc2D wrong'
                    sys.exit()
                r=zeros(n,'f')
                for i in xrange(n):
                    r[i] = self._call_val_bisect(posx[i],posy[i])
                return r
            '''
        else :
            return self._call_val_bisect(posx, posy)
        
        
    def _call_val(self,posx, posy):
        #frist find the index where the data is
        if posx<self.grid[0][0] or posx > self.grid[0][-1]  \
                or posy<self.grid[1][0] or posy > self.grid[1][-1] :
            print 'ERROR: value of func asked in ', posx,',',posy \
                    , 'outside of domain [', self.grid[0][0],',',self.grid[0][-1],']x[' \
                    , self.grid[1][0],',',self.grid[1][-1],']x['
            sys.exit()
        xpos=0; ypos=0
        for gridx in self.grid[0][1:] :
            if posx<=gridx:
                break
            xpos +=1
        for gridy in self.grid[1][1:] :
            if posy<=gridy:
                break
            ypos +=1
        #interpolate
        return self._call_val_data(self.grid[0][xpos],self.grid[0][xpos+1], 
                                self.grid[1][ypos],
                                self.grid[1][ypos+1],
                                posx,posy,
                                self.val[xpos+1,ypos+1],
                                self.val[xpos,ypos+1],
                                self.val[xpos,ypos],
                                self.val[xpos+1,ypos])

    def _call_val_data(self,gridx,gridxp,gridy,gridyp,posx,posy, 
                            val11,val01,val00,val10):
        x0 = (gridx + gridxp) /2.
        y0 = (gridy + gridyp) /2.
        a  = gridxp - gridx
        b  = gridyp - gridy
        xa = 2. * (posx-x0)/a
        ya = 2. * (posy-y0)/b
        phi_1 =  (1+xa)*(1+ya)/4.
        phi_2 =  (1-xa)*(1+ya)/4.
        phi_3 =  (1-xa)*(1-ya)/4.
        phi_4 =  (1+xa)*(1-ya)/4.
        return val11 * phi_1 + val01 * phi_2 + val00 * phi_3 + val10 * phi_4

    def _call_val_bisect(self,posx,posy):
        #first find the index where the data is
        if posx<self.grid[0][0] or posx > self.grid[0][-1] \
                or posy<self.grid[1][0] or posy > self.grid[1][-1] :
            print 'ERROR: value of func asked in ', posx,',',posy \
                    , 'outside of domain [', self.grid[0][0],',',self.grid[0][-1],']x[' \
                    , self.grid[1][0],',',self.grid[1][-1],']x['
            sys.exit()
        xpos = self.lastxpos
        ypos = self.lastypos
        if self.grid[0][xpos] <= posx and self.grid[0][xpos+1]>= posx \
         and self.grid[1][ypos] <= posy and self.grid[1][ypos+1]>= posy:
            return self._call_val_data(self.grid[0][xpos],self.grid[0][xpos+1], 
                            self.grid[1][ypos],self.grid[1][ypos+1],
                            posx,posy, self.val[xpos+1,ypos+1],self.val[xpos,ypos+1],
                            self.val[xpos,ypos],self.val[xpos+1,ypos])
        xpos = bisect_left(self.grid[0],posx)
        xpos -= 1
        if xpos == -1 : xpos = 0
        self.lastxpos = xpos
        ypos = bisect_left(self.grid[1],posy)
        ypos -= 1
        if ypos == -1 : ypos = 0
        self.lastypos = ypos
        #interpolate
        return self._call_val_data(self.grid[0][xpos],self.grid[0][xpos+1],
                            self.grid[1][ypos],self.grid[1][ypos+1],
                            posx,posy, self.val[xpos+1,ypos+1],self.val[xpos,ypos+1],
                            self.val[xpos,ypos],self.val[xpos+1,ypos])

    def _call_val_weave_grid(self,posx,posy):
        # input is an array posx and an array posy, being x grid and y grid
        #returns a matrix with the grid of results xy
        #   No check on data integrity here
        
        #there is a bug somewhere that if I just pass posx, or xvalarray=posx to weave
        # xvalarray is actually empty. I think this is because it references data, instead of being copy
        # it defenitely is a weave bug as it works for posy !
        #A workaround, make a deepcopy:
        xvalarray =  zeros(len(posx),'float')
        xvalarray[:] = posx[:]
        yvalarray =  zeros(len(posy),'float')
        yvalarray[:] = posy[:]
        
        nrx   = len(posx)
        nry   = len(posy)
        result = zeros(nrx*nry, float)
        xvals= self.grid[0]
        nxval= len(self.grid[0])
        xvalmin = self.grid[0][0]
        xvalmax = self.grid[0][-1]
        yvals= self.grid[1]
        nyval= len(self.grid[1])
        yvalmin = self.grid[1][0]
        yvalmax = self.grid[1][-1]
        zvals= self.val  #zval must be numpy array !
        assert(type(xvals)==type([]))
        assert(type(yvals)==type([]))
        # for output, add eg: printf(" %g %g", zval11, result[k*nrx+j]); 
        code = r"""
#line 512 "_call_val_weave_grid"
double xvalsi,yvalsi,xvalsip, yvalsip,xval,yval;
double zval00,zval01,zval10,zval11,x0,y0,a,b,xa,ya,phi_1,phi_2,phi_3,phi_4;
int indexx,indexy, tmpx, tmpy,ind;
int ok;
indexx=0;
indexy=0;
tmpx=0;
tmpy=0;
ind=0;
ok=1;
for (int k=0; k<nrx; k++) {
for (int j=0; j<nry; j++) {
    xval = xvalarray[k];
    yval = yvalarray[j];
    if (xval < xvalmin || xval > xvalmax || yval <yvalmin || yval > yvalmax ) {
        ok=0;
        j=nry;
        k=nrx;
        printf("ERROR GridUtils.py, weave code");
        }
    if (xval>= py_to_float(PyList_GetItem(xvals,tmpx),"xvalt1") && xval<= py_to_float(PyList_GetItem(xvals,tmpx+1),"xvalt2")) {
        indexx=tmpx;
        } else if (tmpx+1 < nxval-1 &&  xval>= py_to_float(PyList_GetItem(xvals,tmpx+1),"xvalt3") && xval<= py_to_float(PyList_GetItem(xvals,tmpx+2),"xvalt4")) {
        indexx=tmpx+1;
        tmpx = indexx;
        } else {
            for (int i=1; i<nxval; i++){
              if (xval <= py_to_float(PyList_GetItem(xvals,i),"xvalsip")) {
                  indexx = i-1;
                  tmpx=indexx;
                  i = nxval;
                  }
              }
        }
    if (yval>= py_to_float(PyList_GetItem(yvals,tmpy),"yvalt1") && yval<= py_to_float(PyList_GetItem(yvals,tmpy+1),"yvalt2")) {
        indexy=tmpy;
        } else if (tmpy+1 < nyval-1 &&  yval>= py_to_float(PyList_GetItem(yvals,tmpy+1),"yvalt3") && yval<= py_to_float(PyList_GetItem(yvals,tmpy+2),"yvalt4")) {
        indexy=tmpy+1;
        tmpy=indexy;
        } else {
            for (int i=1; i<nyval; i++){
              if (yval <= py_to_float(PyList_GetItem(yvals,i),"yvalsip")) {
                  indexy = i-1;
                  tmpy=indexy;
                  i = nyval;
                  }
              }
        }
    
    xvalsip= py_to_float(PyList_GetItem(xvals,indexx+1),"xvalsip");
    xvalsi = py_to_float(PyList_GetItem(xvals,indexx),"xvalsi");
    yvalsi = py_to_float(PyList_GetItem(yvals,indexy),"yvalsi");
    yvalsip = py_to_float(PyList_GetItem(yvals,indexy+1),"yvalsip");
    zval11 = zvals[(indexx+1)*nyval +indexy+1];
    zval01 = zvals[indexx*nyval +indexy+1];
    zval00 = zvals[indexx*nyval +indexy];
    zval10 = zvals[(indexx+1)*nyval+indexy];
    //printf(" %d %d %g %g", indexx, indexy,zval00,zval10); 
    x0 = (xvalsi + xvalsip) /2.;
    y0 = (yvalsi + yvalsip) /2.;
    a  = xvalsip - xvalsi;
    b  = yvalsip - yvalsi;
    xa = 2. * (xval-x0)/a;
    ya = 2. * (yval-y0)/b;
    phi_1 =  (1+xa)*(1+ya)/4.;
    phi_2 =  (1-xa)*(1+ya)/4.;
    phi_3 =  (1-xa)*(1-ya)/4.;
    phi_4 =  (1+xa)*(1-ya)/4.;
    ind=k*nry+j;
    result[ind] = zval11 * phi_1 + zval01 * phi_2 + zval00 * phi_3 + zval10 * phi_4,"res";
}
}
return_val = ok;
"""
        # compiler keyword only needed on windows with MSVC installed
        msg = weave.inline_tools.inline(code,
                           ['nrx', 'nry', 'xvalarray', 'yvalarray', 'nxval', 'xvals', 
                            'nyval', 'yvals', 'zvals', 'result'
                            ,'xvalmin','xvalmax','yvalmin','yvalmax']
                           #,type_converters = converters.blitz
                           ,verbose=1)
        if not msg == 1 :
            print 'ERROR error in _call_val_weave_grid, out of bounds, msg=' ,msg
            print 'grid x', posx
            print 'grid y', posy
            print 'grid ', self.grid
            sys.exit()
        #we make result a two-dim array:
        result = resize(result,(nrx,nry))
        return result
        
    def _call_val_weave(self,posx,posy):
        # input is an numpy array posx of same length as posy
        #   No check on data integrity here
        
        #there is a bug somewhere that if I just pass posx, or xvalarray=posx to weave
        # xvalarray is actually empty. I think this is because it references data, instead of being copy
        # it defenitely is a weave bug as it works for posy !
        #A workaround, make a deepcopy:
        xvalarray =  zeros(len(posx),'float')
        xvalarray[:] = posx[:]
        yvalarray =  zeros(len(posy),'float')
        yvalarray[:] = posy[:]
        
        nr   = len(posx)
        result = zeros(nr, float)
        xvals= self.grid[0]
        nxval= len(self.grid[0])
        xvalmin = self.grid[0][0]
        xvalmax = self.grid[0][-1] 
        yvals= self.grid[1]
        nyval= len(self.grid[1])
        yvalmin = self.grid[1][0]
        yvalmax = self.grid[1][-1]
        zvals= self.val  #zval must be numpy array !
        assert(type(xvals)==type([]))
        assert(type(yvals)==type([]))
        code = r"""
#line 601 "_call_val_weave"
double xvalsi,yvalsi,xvalsip, yvalsip,xval,yval;
double zval00,zval01,zval10,zval11,x0,y0,a,b,xa,ya,phi_1,phi_2,phi_3,phi_4;
int indexx,indexy, tmpx, tmpy;
int ok;
indexx=0;
indexy=0;
tmpx=0;
tmpy=0;
ok=1;
for (int j=0; j<nr; j++) {
    xval = xvalarray[j];
    yval = yvalarray[j];
    if (xval < xvalmin || xval > xvalmax || yval <yvalmin || yval > yvalmax ) {
        ok=0;
        j=nr;
        for (int k=0; k<nr;k++) {
            result[k]=0.00001;
            }
        } else 
    {
    if (xval>= py_to_float(PyList_GetItem(xvals,tmpx),"xvalt1") && xval<= py_to_float(PyList_GetItem(xvals,tmpx+1),"xvalt2")) {
        indexx=tmpx;
        } else if (tmpx+1 < nxval-1 &&  xval>= py_to_float(PyList_GetItem(xvals,tmpx+1),"xvalt3") && xval<= py_to_float(PyList_GetItem(xvals,tmpx+2),"xvalt4")) {
        indexx=tmpx+1;
        tmpx = indexx;
        } else {
            for (int i=1; i<nxval; i++){
              if (xval <= py_to_float(PyList_GetItem(xvals,i),"xvalsip")) {
                  indexx = i-1;
                  tmpx=indexx;
                  i = nxval;
                  }
              }
        }
    if (yval>= py_to_float(PyList_GetItem(yvals,tmpy),"yvalt1") && yval<= py_to_float(PyList_GetItem(yvals,tmpy+1),"yvalt2")) {
        indexy=tmpy;
        } else if (tmpy+1 < nyval-1 &&  yval>= py_to_float(PyList_GetItem(yvals,tmpy+1),"yvalt3") && yval<= py_to_float(PyList_GetItem(yvals,tmpy+2),"yvalt4")) {
        indexy=tmpy+1;
        tmpy=indexy;
        } else {
            for (int i=1; i<nyval; i++){
              if (yval <= py_to_float(PyList_GetItem(yvals,i),"yvalsip")) {
                  indexy = i-1;
                  tmpy=indexy;
                  i = nyval;
                  }
              }
        }
    
    xvalsip= py_to_float(PyList_GetItem(xvals,indexx+1),"xvalsip");
    xvalsi = py_to_float(PyList_GetItem(xvals,indexx),"xvalsi");
    yvalsi = py_to_float(PyList_GetItem(yvals,indexy),"yvalsi");
    yvalsip = py_to_float(PyList_GetItem(yvals,indexy+1),"yvalsip");
    zval11 = zvals[(indexx+1)*nyval +indexy+1];
    zval01 = zvals[indexx*nyval +indexy+1];
    zval00 = zvals[indexx*nyval +indexy];
    zval10 = zvals[(indexx+1)*nyval+indexy];
    x0 = (xvalsi + xvalsip) /2.;
    y0 = (yvalsi + yvalsip) /2.;
    a  = xvalsip - xvalsi;
    b  = yvalsip - yvalsi;
    xa = 2. * (xval-x0)/a;
    ya = 2. * (yval-y0)/b;
    phi_1 =  (1+xa)*(1+ya)/4.;
    phi_2 =  (1-xa)*(1+ya)/4.;
    phi_3 =  (1-xa)*(1-ya)/4.;
    phi_4 =  (1+xa)*(1-ya)/4.;
    result[j] = zval11 * phi_1 + zval01 * phi_2 + zval00 * phi_3 + zval10 * phi_4;
    }
}
return_val = ok;
"""
        # compiler keyword only needed on windows with MSVC installed
        msg = weave.inline_tools.inline(code,
                           ['nr', 'xvalarray', 'yvalarray', 'nxval', 'xvals', 
                            'nyval', 'yvals', 'zvals', 'result'
                            ,'xvalmin','xvalmax','yvalmin','yvalmax']
                           #,type_converters = converters.blitz
                           ,verbose=1)
        if not msg == 1 :
            print 'ERROR error in _call_val_weave, out of bounds' 
            print 'x bounds: ',xvalmin,xvalmax,
            print 'x data:', xvalarray
            print 'y bounds: ',yvalmin,yvalmax,
            print 'y data:', yvalarray
            print ' In case of optimization: reduce stepsize.'
            raise DomainError, [] 
            sys.exit()
        return result
                             
class OptPath:
    ''' Optimization path. 
        A lists of points defines a path on a 2D grid. A left and right size 
        define the zone in which to optimize
    '''
    def __init__(self, path, left, right, down, up):
        ''' Give a path as [[0,0],[1,1], [2,1]]
            and double values of left, right, up, down
            The path must be monotone in x (x[i] <= x[i+1] for all x or >=)
             AND monotone in y. Split paths in parts that are monotone if this is not satisfied!!
        '''
        xval = [x[0] for x in path]
        yval = [x[1] for x in path]
        #check monotone 
        signx = self._test_monotone(xval)
        if not signx :
            raise ValueError, 'Data is not monotone in x'
        signy = self._test_monotone(yval)
        if not signy :
            raise ValueError, 'Data is not monotone in y' 
        if signx =='<=' :
            self.leftright = GridFunc1D([xval], yval)
            self.minx = xval[0]
            self.maxx = xval[-1]
        else :
            self.leftright = GridFunc1D([xval.reverse()], yval.reverse())
            self.minx = xval[-1]
            self.maxx = xval[0]
        if signy =='<=':
            self.downup    = GridFunc1D([yval], xval)
            self.miny = yval[0]
            self.maxy = yval[-1]
        else :
            self.downup    = GridFunc1D([yval.reverse()], xval.reverse())
            self.miny = yval[-1]
            self.maxy = yval[0]
                    
        self.left=left
        self.right=right
        self.up=up
        self.down=down
        
    def _test_monotone(self,xval):
        #check monotone 
        signx = None
        xarr = N.array(xval)
        #xarr and t must be completely equal type !
        t    = zeros(len(xval),type(xarr[0]))
        t[:-1]=xval[1:]
        t[-1] =xval[-1]
        testsign = (xarr <= t)
        signx = '<='
        for test in testsign :
            if not test :
                signx = None
                break
        if signx == None :
            testsign = (xarr >= t)
            signx = '>='
            for test in testsign :
                if not test :
                    signx = None
                    break
        return signx
        
    def near(self, point):
        ''' Determine if a point is near the optimization path
            point = [double, double]
            returns true or false
        '''
        if point[0] < self.minx - self.left or point[0]>self.maxx + self.right \
                 or point[1] < self.miny - self.down or point[1]>self.maxy + self.up:
            return False
        if point[0] < self.minx :
            valud = self.leftright(self.minx)
        elif point[0] > self.maxx :
            valud = self.leftright(self.maxx)
        else :
            valud = self.leftright(point[0])
        if valud + self.up >= point[1] and valud - self.down <= point[1] :
                return True
        if point[1] < self.miny :
            vallr = self.downup(self.miny)
        elif point[1]>self.maxy :
            vallr = self.downup(self.maxy)
        else :
            vallr = self.downup(point[1])
        if vallr + self.right >= point[0] and vallr - self.left <=point[0]:
            return True
        return False
    
    def optimize(self,point):
        ''' is it necessary to optimize in a point? 
            Yes if the point is near, false otherwise
        '''
        return self.near(point)
    
    def plotongrid(self,xvals,yvals, asciart = True):
        ''' Given a grid: [xval[i],yval[j]], go over all entries in the grid and print to
            output where we optimize
        '''
        if not isinstance(yvals,list) :
            yvalsrev = yvals.tolist()
            yvalsrev.reverse()
        else:
            yvalsrev = copy.deepcopy(yvals)
            yvalsrev.reverse()

        for yval in yvalsrev :
            line = ''
            for xval in xvals :
                if self.near([xval,yval]) :
                    line += 'x '
                else :
                    line += '  '
            print line
            
                
        
                                    
def _test_gridfunc(): 
    
    def plot(uL, uR, step, func):
        import pylab as pl
        u = pl.arange(uL, uR+step/2., step)
        fu = func(u)
        #print u
        #print Du
        pl.plot(u, fu)
        pl.xlabel('x val')
        pl.ylabel('func val')
        #pl.title('Diffusion')
        pl.grid(True)
        pl.show()
    def plot2D(uL, uR, stepu, vL, vR, stepv, func):
        ut = N.arange(uL, uR+stepu/2., stepu)
        vt = N.arange(vL, vR+stepv/2., stepv)
        fu = func(ut,vt)
        print ut,vt
        print fu
        fu = resize(fu,(size(ut),size(vt)))
        #plot with gnuplot
        import Gnuplot
        g  = Gnuplot.Gnuplot(persist=1)
        g.title('test1')
        g.splot(Gnuplot.GridData(fu ,ut, vt,binary=True))
        
    print 'test in 1D '
    grid = N.arange(0,3.01,0.5,float)
    grid= [grid.tolist(),]
    grid_small = N.arange(0,3.001,0.01,float)
    grid_small= [grid_small.tolist(),]
    val1 = N.exp(grid[0]).tolist()
    val2 = N.cos(grid[0]).tolist()
    a=GridFunc1D(grid,val1)
    b=GridFunc1D(grid,val2)
    val1_small = N.exp(grid_small[0]).tolist()
    val2_small = N.cos(grid_small[0]).tolist() 
    a_small=GridFunc1D(grid_small,val1_small)
    b_small=GridFunc1D(grid_small,val2_small)
    #plota=[]
    #plotb=[]
    #for i in range(100) :
    #    plota += [a(i*3./100.)]
    #    plotb += [b(i*3./100.)]
    #print plota
    #print plotb
    plot(0,3,0.02,a)
    plot(0,3,0.02,a_small)
    plot(0,3,0.02,b)
    plot(0,3,0.02,b_small)
    
    print 'test in 2D '
    ut = N.mgrid[0:3:7j].tolist()
    vt = N.mgrid[0:3:7j].tolist()
    fu = []
    for i in range(size(ut)):
        for j in range(size(vt)):
            fu.append( sin(ut[i]*vt[j]))
    fu = resize(fu,(size(ut),size(vt)))
    print fu   
    a=GridFunc2D([ut,vt],N.asarray(fu))
    #plota=[]
    #plotb=[]
    #for i in range(100) :
    #    plota += [a(i*3./100.)]
    #    plotb += [b(i*3./100.)]
    #print plota
    #print plotb
    plot2D(0,3,0.2,0,3,0.2,a)
    ut_small = N.mgrid[0:3:100j].tolist()
    vt_small = N.mgrid[0:3:100j].tolist()
    fu_small = []
    for i in range(size(ut_small)):
        for j in range(size(vt_small)):
            fu_small.append( sin(ut_small[i]*vt_small[j]))
    fu_small = resize(fu_small,(size(ut_small),size(vt_small)))   
    a_small=GridFunc2D([ut_small,vt_small],N.asarray(fu_small))
    plot2D(0,3,0.2,0,3,0.2,a_small)
    print ' test of the three calling methods : '
    print 'GridFunc2D test point    ' , a_small(1.5,1.5), type(a_small(1.5,1.5))
    print 'GridFunc2D list-gridout  ' , a_small([1.,1.5,2.],[1.5,1.8]), type(a_small([1.,1.5,2],[1.5,1.8])) , type(a_small([1,1.5,2],[1.5,1.8])[1])
    print 'GridFunc2D list-nogrid   ' , a_small([1.,1.5,2.],[1.5,1.8,1.8],False), type(a_small([1.,1.5],[1.5,1.8],False)), type(a_small([1,1.5],[1.5,1.8],False)[1])
    
def _test_optpath():
    a = OptPath([[0,0],[1,1],[2,1]],0.1,0.2,0.1,0.2)
    print 'True?',a.near([-0.1,0])
    print 'True?',a.near([0.2,0])
    print 'True?',a.near([0,0])
    print 'False?',a.near([0.21,0])
    print 'True?',a.near([0.5,0.6])
    print 'True?',a.near([1.5,0.95])
    print 'True?',a.near([1.5,1.15])
    print 'False?',a.near([1.5,1.20001])
    
    xgrid = N.mgrid[0.:2.:21j].tolist()
    ygrid = N.mgrid[-0.5:1.5:21j].tolist()
    a.plotongrid(xgrid, ygrid)
    b = OptPath([[0,0],[1,1],[2,1]],0.2,0.2,0.2,0.2)
    b.plotongrid(xgrid, ygrid)


if __name__ == "__main__":
    _test_gridfunc()
    _test_optpath()                                    
                            