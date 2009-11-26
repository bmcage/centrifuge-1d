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

"""
    Read data from command line arguments
    Open files with experiment data
    Construct a mesh
    Put data in a mesh
"""
#-------------------------------------------------------------------------
#
# Python modules
#
#-------------------------------------------------------------------------
import sys, os, shutil
import numpy as N
import scipy as S
from optparse import OptionParser
import Gnuplot
from math import *

#-------------------------------------------------------------------------
#
# local modules
#
#-------------------------------------------------------------------------
import GridUtils
import UtilsBm
import centrifuge
import satunsat

#-------------------------------------------------------------------------
#
# Constants
#
#-------------------------------------------------------------------------

eps = 1e-6

#-------------------------------------------------------------------------
#
# the main prog
#
#-------------------------------------------------------------------------

def main():
    
    parser = OptionParser()
    add_parser_options(parser)
    
    options, args = parser.parse_args(sys.argv[1:])
    
    # make mesh on [0,1]
    grid = GridUtils.create_grid(options.nrp, 0., 1., "borderrefined", 
                                 options.refineval)
    #read in initial conc
    startgrid = options.leftposL
    endgrid = options.L
    timeinit = 0.;
    meshinitpos = []               
    pinit  = []
    initfunction = eval(options.init)
    for pos in grid : 
        pp = pos * (endgrid-startgrid) + startgrid
        meshinitpos.append(pp)  
        pinit.append(initfunction(pp)) # * centrifuge.rho_f * centrifuge.g )

    
    # make data for the values of initial cond over the grid
    measval = S.zeros(options.nrp , float)
   
    # projec measurements onto grid
    measval[:] = pinit[:]
                        
    dir = options.dir
    if os.path.isdir(dir):
        shutil.rmtree(dir)
    os.mkdir(dir)
    
    # we plot with gnuplot the initial condition and its projection to a grid 
    plot_init_cond(os.path.join(dir, 'init.ps'), meshinitpos, pinit)
    

    print '  ---- Starting direct solution ----'
    #we have a centrifuge flow problem, over the grid, with init conc measval
    # and extra data of parameters in the options class.
    data = centrifuge.DataCentrifuge(options)
    print data.data()
    if data.inittype() =="satinflow":
        Model = centrifuge.CentrifugeFlow(grid, measval, data, dir=options.dir)
    elif data.inittype() =="satunsatinflow":
        Model = satunsat.SatUnsatFlow(grid, measval, data, dir=options.dir)
    else:
        assert False, data.inittype()+ ' not supported'
    Model.solve()
    
    #print the mass out
    Model.showmass(Model.get_outputdir())
    
    #plot solution
    if options.saveoutput :
        #try :
            Model.save_output(Model.get_outputdir())
        #except AttributeError:
        #    pass
    Model.plot_sol(Model.get_outputdir())
    
    #print out residual
    res = Model.residuals()
    if res is not None:
        print 'Residual of the Model compared to experiment: ',res
    
def add_parser_options(parser):
    """
    default options, and add the options anneal programm recognises
    """
    #initial concentration defaults
    parser.set_defaults(init=  "lambda x: 0")
    #grid defaults
    parser.set_defaults(nrp=101, refineval=1.)
    #problem defaults
    parser.set_defaults(L=0.06)
    parser.set_defaults(r=0.1)
    parser.set_defaults(g="lambda t: 1.")
    parser.set_defaults(omega="lambda t: 10.")
    parser.set_defaults(leftpos=0.5)
    parser.set_defaults(leftposL=0.6)
    parser.set_defaults(th_p1=0.01)
    parser.set_defaults(th_p2=0.01)
    parser.set_defaults(k_p=1e-13)
    parser.set_defaults(n_p=0.5)
    parser.set_defaults(temp=20.)
    parser.set_defaults(k=1.2e-19)
    parser.set_defaults(mu=1.003e-3)
    parser.set_defaults(rho_f=998.2071)
    parser.set_defaults(grav=9.80665)
    parser.set_defaults(n0=0.3)
    parser.set_defaults(alpha=7.86e-9)
    parser.set_defaults(alphagen=-0.0189)
    parser.set_defaults(ngen=2.81)
    parser.set_defaults(theta_r=0.02)
    parser.set_defaults(theta_s=0.4)
    parser.set_defaults(alphagen_p=-0.0189)
    parser.set_defaults(ngen_p=2.81)
    parser.set_defaults(theta_r_p=0.02)
    parser.set_defaults(theta_s_p=0.4)
    parser.set_defaults(eps=0.0001)
    parser.set_defaults(BC="dirhomdir")
    parser.set_defaults(inittype="satinflow")
    parser.set_defaults(s0=0.001)
    parser.set_defaults(pcalc="p")
    parser.set_defaults(endTime = 1600.0, nrsteps = 1, outputtime=200)
    parser.set_defaults(dir="test")
    
    parser.add_option('', '--saveoutput', action="store_true", 
                      dest="saveoutput", default=False,
                      help="Save the output as obtained over the outputtime "
                           "timesteps")
    parser.add_option('', '--init', dest='init',
                      help="initial concentration is a function in space")
    parser.add_option('-n', '--numberGridPoints', type='int', dest='nrp',
                      help='number of gridpoints' )
    parser.add_option('','--refineval',type='float',dest='refineval',
                      help='how many times at border more refined than in '
                           'center of grid' )
    parser.add_option('-L','--length',type='float',dest='L',
                      help='Length of the sample' )
    parser.add_option('-g', '--gravity', dest='g', help='earth gravity as a '
                            'function in t')
    parser.add_option('', '--omega', dest='omega', help='rotational speed of '
                            'centrifuge in function of t')
    parser.add_option('', '--leftpos',type='float', dest='leftpos', help='position of '
                            'begin open container')
    parser.add_option('-r', '--leftposL',type='float', dest='leftposL', help='position of '
                            'begin sample')
    parser.add_option('', '--th_p1',type='float', dest='th_p1', help='thickness porous '
                            'stone to left (side axis)')
    parser.add_option('', '--th_p2',type='float', dest='th_p2', help='thickness porous '
                            'stone to right ')
    parser.add_option('', '--k_p',type='float', dest='k_p', help='intrinsic permeability '
                        'porous stone')
    parser.add_option('', '--n_p',type='float', dest='n_p', help='porosity '
                        'porous stone')
    parser.add_option('', '--temp',type='float', dest='temp', help='temperature of sample')
    parser.add_option('-k', '--permea',type='float', dest='k', help='intrinsic permeability'
                            ' in sample')
    parser.add_option('', '--mu',type='float', dest='mu', 
                        help='viscosity of water at temp')
    parser.add_option('', '--rho_f',type='float', dest='rho_f',
                        help='density of water at temp')
    parser.add_option('', '--grav',type='float', dest='grav', help='gravitational '
                            'accelleration g, default=9.80665 m/s^2')
    parser.add_option('', '--n0', type='float', dest='n0', 
                            help='reference porosity in sample')
    parser.add_option('', '--alpha', type='float', dest='alpha', 
                            help='compressibility porous medium')
    parser.add_option('', '--alphagen', type='float', dest='alphagen', 
                            help='alpha parameter in Van Genuchten')
    parser.add_option('', '--ngen', type='float', dest='ngen', 
                            help='n parameter in Van Genuchten')
    parser.add_option('', '--theta_r', type='float', dest='theta_r', 
                            help='theta_r parameter in Van Genuchten')
    parser.add_option('', '--theta_s', type='float', dest='theta_s', 
                            help='theta_s parameter in Van Genuchten')
    parser.add_option('', '--alphagen_p', type='float', dest='alphagen_p', 
                            help='alpha parameter in Van Genuchten for porous stone')
    parser.add_option('', '--ngen_p', type='float', dest='ngen_p', 
                            help='n parameter in Van Genuchten for porous stone')
    parser.add_option('', '--theta_r_p', type='float', dest='theta_r_p', 
                            help='theta_r parameter in Van Genuchten for porous stone')
    parser.add_option('', '--theta_s_p', type='float', dest='theta_s_p', 
                            help='theta_s parameter in Van Genuchten for porous stone')
    parser.add_option('', '--eps', type='float', dest='eps', 
                            help='eps parameter modelling the storativity in the model')
    parser.add_option('', '--expfile', action='append', dest='expfile',
                      help='filename where experiment conc are found, several '
                           'experiments can be given')    
    parser.add_option('-T', '--endTime', type='float', dest='endTime',
                      help='end time of the model' )
    parser.add_option('-t', '--timesteps', type='int', dest='nrsteps',
                      help='number of steps over time interval, needed in '
                           'the ODE solver' )
    parser.add_option('', '--outputtime', type='int', dest='outputtime',
                      help='time at wich output of pic of solution must '
                           'happen' )
    #parser.add_option('-K','--hydcond',type='float',dest='K',
    #                  help='hydraulic conductivity in sample' )
    parser.add_option('', '--pcalc', dest='pcalc', help='input given as '
                'pressures p1,p2, or given in height of water table: p or h')
    parser.add_option('', '--BC', dest='BC',
                      help='Type of boundary  possible: '
                           '"dirhomdir", ...') 
    parser.add_option('', '--inittype', dest='inittype',
                      help='Type of initial values at start, possible: '
                           '"satinflow" : saturated, water going in (default)'
                           '"satunsatinflow", ...') 
    parser.add_option('','--p1',type='float',dest='p1',
                      help='pressure left BC (at 0)' )
    parser.add_option('','--p2',type='float',dest='p2',
                      help='pressure right BC (at L)' )
    parser.add_option('','--s0',type='float',dest='s0',
                      help='initial water table position relative to '
                            'begin container' )
    parser.add_option('','--h2',type='float',dest='h2',
                      help='head right BC (at L)' )
    parser.add_option('', '--dir', dest='dir',
                      help='The directory where output files will be stored')

def plot_init_cond(filename, mesh, pressure):
    g = Gnuplot.Gnuplot(persist=0)
    pSi     = Gnuplot.Data(mesh, pressure, with_='points', 
                           title='Initial condition')
    g.title('Initial pressure distribution')
    g.xlabel('x (m)')
    g.ylabel('p (Pa)')
    g.plot(pSi)
    #make postscript file with plot
    g.hardcopy(filename=filename, enhanced=1,
               mode='eps', color=0, fontname='Times-Roman', fontsize=28)
            
    g = Gnuplot.Gnuplot(persist=0)

if __name__ == "__main__":
    main()

