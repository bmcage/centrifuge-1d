from __future__ import division, print_function

from ..shared.functions import lagrangian_derivative_coefs
from numpy import linspace, power, empty
from ..shared.consolidation import (create_CON, CON_SLURRY, CON_GOMPERTZ,
                                    CON_FREEFORM)

def dtype_deps(cfg):
    dtype = cfg.get_value('dtype')
    result = []
    if dtype == 1: pass
    elif dtype in [2,3]: result = ['k_dx']

    return result

PARENTAL_MODULES = ['base']

CONFIG_OPTIONS = ['inner_points', 'dtype',
                  ('con_type', CON_SLURRY),
                  ('con_max_refine', 0),
                  (lambda cfg: cfg.get_value('con_type') == CON_SLURRY,
                   ['a', 'b', 'c', 'd']),
                  (lambda cfg: cfg.get_value('con_type') in [CON_GOMPERTZ,],
                   ['a', 'b', 'c', 'd', 'cc']),
                  (lambda cfg: cfg.get_value('con_type') in [CON_FREEFORM],
                   [('ei', None), ('si', None), ('ki', None), ('eiadd', None)]),
                  'porosity',
                  'estimate_zp0',
                  ('L_atol', 1e-8),
                  dtype_deps,
                  # dependent
                  (lambda cfg: cfg.get_value('fl1') > 0.0,
                      ['fp1'], [('fp1', -1.0)]),
                  (lambda cfg: cfg.get_value('fl2') > 0.0,
                      ['fp2'], [('fp2', -1.0)]),
                  #
                  'rb_type',
                  # dependent
                  (lambda cfg: cfg.get_value('rb_type') == 2,
                    ['h_last']),
                  (lambda cfg: cfg.get_value('rb_type') == 3,
                    ['dip_height']),
                  'h_last',
                  'l0',
                  'wl0',
                  'density_s', #density sample in g/(cm^3)
                  ('numerfact_e0', 0.999)
                 ]

INTERNAL_OPTIONS = ['m', 'y', 'y12', 'dy', 'alpha', 'ldc1', 'ldc2', 'ldc3',
                    'k_dx', 'wm0', 'CON',
                    'first_idx', 'last_idx', 'wl_idx', 'L_idx',
                    'mass_in_idx', 'mass_out_idx',
                    'z_size', 'gamma_w', 'gamma_s', 'e0']

EXCLUDE_FROM_MODEL = ['dtype']

PROVIDE_OPTIONS = []

OPTIONS_ITERABLE_LISTS = ['porosity']


def adjust_cfg(cfg):
    #specific weight water in g/(s cm^2)
    cfg.set_value('gamma_w', cfg.get_value('density')*cfg.get_value('g'))
    #specific weight sample in g/(s cm^2)
    cfg.set_value('gamma_s', cfg.get_value('density_s')*cfg.get_value('g'))

    # Discretization
    inner_points = cfg.get_value('inner_points')

    discretization_type = cfg.get_value('dtype')
    if discretization_type == 1:   # linear discretization
        y = linspace(0, 1, inner_points + 2)
    elif discretization_type in [2,3]: # L= a+ka+(k^2)a+...+(k^inner_points)a
        # L=1 (as we use transformed interval <0,1>)
        # L = a*[(1-k^(inner_points +1))/(1-k)]
        k = cfg.get_value('k_dx')
        a=(1-k)/(1-power(k, inner_points+1))
        y= empty([inner_points+2, ])
        y[0] = 0.0; y[-1] = 1.0
        for i in range(1, inner_points+1):
            y[i] = y[i-1] + a
            a = a*k
        if discretization_type == 3:
            # invert it
            tmp = y[::-1]
            y[:] = 1. - tmp[:]
    else:
        print('Unsupported discretization type:', discretization_type)
        exit(1)

    #porosity and void ratio
    por = cfg.get_value('porosity')
    if not (0<por<1):
        print ('Porosity must be a value between 0 and 1. Given:', por)
        exit(1)
    e0 = por/(1-por)
    cfg.set_value('e0', e0)
    print ('Consolidation: Calculated initial void ratio is', cfg.get_value('e0'))
    ksguess = cfg.get_value('ks')
    ks = ksguess
    if cfg.get_value('con_type') == CON_SLURRY:
        ks = (1+e0)*(cfg.get_value('c')+cfg.get_value('d')*e0)
        cfg.set_value('ks', ks)
    print ('Consolidation: Your guessed ks', ksguess, 'has been changed into calculated', ks)

    # Determine consolidation curve model used, all data is now available
    cfg.set_value('CON', create_CON(cfg))

    cfg.set_value('y', y)
    cfg.set_value('y12', (y[1:]+y[:-1])/2.)

    dy = y[1:]-y[:-1]
    alpha = empty([len(dy)+1, ])
    alpha[0] = 0.
    alpha[1:] = dy
    cfg.set_value('dy', dy)
    cfg.set_value('alpha', alpha)

    ldc1, ldc2, ldc3 = lagrangian_derivative_coefs(dy)
    cfg.set_value('ldc1', ldc1)
    cfg.set_value('ldc2', ldc2)
    cfg.set_value('ldc3', ldc3)

    inner_points = cfg.get_value('inner_points')
    cfg.set_value('sc_max_refine', 0)
    cfg.set_value('first_idx',    0)
    cfg.set_value('last_idx',     inner_points+1)
    cfg.set_value('mass_in_idx',  inner_points+2)
    cfg.set_value('wl_idx',       inner_points+3)
    cfg.set_value('L_idx',        inner_points+4)
    cfg.set_value('mass_out_idx', inner_points+5)
    # total length of 'z' array (discretization points + s1,s2,mass_in,...)
    cfg.set_value('z_size',       inner_points+6)


def check_cfg(cfg):
    if not (not cfg.get_value('wl0') is None or not cfg.get_value('ww0') is None):
        print("One of 'wl0' or 'ww0' parameters must be specified.")
        return False

    return True
