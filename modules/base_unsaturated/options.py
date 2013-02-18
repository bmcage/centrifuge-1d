from __future__ import division

import modules.shared.saturation_curve as mSC
from modules.shared.functions import lagrangean_derivative_coefs
from numpy import linspace, power, empty

def dtype_deps(cfg):
    dtype = cfg.get_value('dtype')
    result = []
    if dtype == 1: pass
    elif dtype == 2: result = ['k_dx']

    return result

PARENTAL_MODULES = ['base']

CONFIG_OPTIONS = ['inner_points', 'dtype', 'n', 'gamma',
                  'porosity',
                  'estimate_zp0',
                  dtype_deps,
                  # dependent
                  (lambda cfg: cfg.get_value('fl1') > 0.0,
                      ['fp1'], [('fp1', -1.0)]),
                  (lambda cfg: cfg.get_value('fl2') > 0.0,
                      ['fp2'], [('fp2', -1.0)])
                 ]

INTERNAL_OPTIONS = ['m', 'y', 'y12', 'dy', 'ldc1', 'ldc2', 'ldc3',
                    'k_dx', 'wm0', 'calc_gc', 'calc_rm', 'calc_wm',
                    'SC']

EXCLUDE_FROM_MODEL = ['dtype']

PROVIDE_OPTIONS = []

OPTIONS_ITERABLE_LISTS = ['porosity']

def adjust_cfg(cfg):
    # Handle depending variables
    SC = mSC.SC_vanGenuchten(cfg.get_value('n'), cfg.get_value('gamma'))
    cfg.set_value('SC', SC)

    # Discretization
    inner_points = cfg.get_value('inner_points')

    discretization_type = cfg.get_value('dtype')
    if discretization_type == 1:   # linear discretization
        y = linspace(0, 1, inner_points + 2)
    elif discretization_type == 2: # L= a+ka+(k^2)a+...+(k^inner_points)a
        # L=1 (as we use transformed interval <0,1>)
        # L = a*[(1-k^(inner_points +1))/(1-k)]
        k = cfg.get_value('k_dx')
        a=(1-k)/(1-power(k, inner_points+1))
        y= empty([inner_points+2, ])
        y[0] = 0.0; y[-1] = 1.0
        for i in range(1, inner_points+1):
            y[i] = y[i-1] + a
            a = a*k
    else:
        print('Unsupported discretization type:', dtype)
        exit(1)

    cfg.set_value('y', y)
    cfg.set_value('y12', (y[1:]+y[:-1])/2.)

    dy = y[1:]-y[:-1]
    cfg.set_value('dy', dy)

    ldc1, ldc2, ldc3 = lagrangean_derivative_coefs(dy)
    cfg.set_value('ldc1', ldc1)
    cfg.set_value('ldc2', ldc2)
    cfg.set_value('ldc3', ldc3)

    cfg.set_value('calc_wm', True)
    cfg.set_value('calc_gc', True)
    cfg.set_value('calc_rm', True)
