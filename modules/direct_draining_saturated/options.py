from __future__ import division

from modules.shared.functions import lagrangean_derivative_coefs
from numpy import linspace, power, asarray, empty

def dtype_deps(cfg):
    dtype = cfg.get_value('dtype')
    result = []
    if dtype == 1: pass
    elif dtype == 2: result = ['k_dx']

    return result

PARENTAL_MODULES = ['base']

CONFIG_OPTIONS = ['inner_points', 'dtype', 'n', 'gamma',
                  'h_init', 'porosity',
                  'calc_gc', 'calc_rm', 'calc_wm',
                  'rb_type',
                  'estimate_zp0',
                  dtype_deps,
                  # dependent
                  (lambda cfg: cfg.get_value('rb_type') == 2,
                    ['h_last']),
                  (lambda cfg: cfg.get_value('rb_type') == 3,
                    ['dip_height', 'h_last',
                     ('s2_0', None), ('s2_atol', 1e-8)]),
                  'h_last',
                  (lambda cfg: cfg.get_value('fl1') > 0.0,
                      ['fp1'], [('fp1', -1.0)]),
                  (lambda cfg: cfg.get_value('fl2') > 0.0,
                      ['fp2'], [('fp2', -1.0)])
                 ]

INTERNAL_OPTIONS = ['m', 'y', 'y12', 'dy', 'ldc1', 'ldc2', 'ldc3',
                    'first_idx', 'last_idx', 's1_idx', 's2_idx',
                    'mass_in_idx', 'mass_out_idx', 'pq_idx', 'z_size',
                    'wm0']

EXCLUDE_FROM_MODEL = ['dtype']

PROVIDE_OPTIONS = []

OPTIONS_ITERABLE_LISTS = ['porosity']

def adjust_cfg(cfg):
    # Handle depending variables
    value = cfg.get_value('n')
    if type(value) == list:
        m = [1.-1./n for n in value]
    else:
        m = 1. - 1./value
    cfg.set_value('m', m)

    # Set array indexes
    inner_points = cfg.get_value('inner_points')

    cfg.set_value('first_idx',    0)
    cfg.set_value('last_idx',     inner_points+1)
    cfg.set_value('mass_in_idx',  inner_points+2)
    cfg.set_value('s1_idx',       inner_points+3)
    cfg.set_value('s2_idx',       inner_points+4)
    cfg.set_value('mass_out_idx', inner_points+5)
    cfg.set_value('pq_idx',       inner_points+6)
    # total length of 'z' array (discretization points + s1,s2,mass_in,...)
    cfg.set_value('z_size',       inner_points+7)

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

    cfg.set_value('calc_wm', cfg.get_value('calc_wm')
                             or cfg.get_value('calc_gc'))
