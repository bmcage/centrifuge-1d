from __future__ import division, print_function
import numpy as np
from collections import OrderedDict

from ..shared.saturation_curve import (SC_vG, SC_FF, SC_FF_BS)

PARENTAL_MODULES = []

CONFIG_OPTIONS = ['inv_init_params', ('optimfn', 'leastsq'),
                  ('transform_params', True), ('untransformed_cov', False),
                  (lambda cfg: cfg.get_value('optimfn') == 'leastsq',
                    ['epsfcn', 'factor',
                     ('xtol', 1.49012e-8), ('ftol', 1.49012e-8)]),
                  (lambda cfg: cfg.get_value('optimfn') == 'fmin_slsqp',
                    ['epsfcn', 'factor',
                     ('xtol', 1.49012e-8), ('max_inv_iter', None),
                     ('disp_inv_conv', True)]),
                  (lambda cfg: (cfg.get_value('optimfn')
                                in ['fmin', 'fmin_powell']),
                    [('xtol', 1e-4), ('ftol', 1e-4), ('max_fev', None),
                     ('max_inv_iter', None), ('disp_inv_conv', True)]),
                  (lambda cfg:
                       cfg.get_value('optimfn') in ['fmin_cg', 'fmin_bfgs'],
                    [('gtol', 1e-5), ('max_inv_iter', None),
                     ('disp_inv_conv', True)])
                  ]

INTERNAL_OPTIONS = ['_transform', '_untransform', 'init_values',
                    '_lbounds', '_ubounds', '_conditions']

EXCLUDE_FROM_MODEL = ['inv_init_params']

def provide_inv_init_params(cfg):
    inv_init_params = cfg.get_value('inv_init_params')

    return list(inv_init_params.keys())

PROVIDE_OPTIONS = [provide_inv_init_params]

OPTIONS_ITERABLE_LISTS = []

def check_cfg(cfg):
    inv_init_params = cfg.get_value('inv_init_params')
    if (not inv_init_params or not type(inv_init_params) is dict):
        print("Option 'inv_init_params' has to be of type dict.")
        return False

    return True

def prior_adjust_cfg(cfg):
    # Assign data from 'inv_init_params' to cfg
    for (name, value) in cfg.get_value('inv_init_params').items():
        if not cfg.get_value(name, not_found='NoTfOuNd') == 'NoTfOuNd':
            print("Parameter '{}' is present in configuration file. It will"
                  " be shadowed with initial value from 'inv_init_params'."
                  .format(name))
        cfg.set_value(name, value[0]) # is a list of (value, range, opts)

def adjust_cfg(cfg):
    # Process 'inv_init_params'
    init_values = OrderedDict()
    lbounds     = {}
    ubounds     = {}
    conditions  = {}

    for (name, value) in cfg.get_value('inv_init_params').items():
        if value is None: continue

        cond = ''
        if len(value) == 1:
            init_value = value
            lbound = -np.inf
            ubound = np.inf
        elif len(value) == 2:
            (init_value, (lbound, ubound)) = value
        elif len(value) == 3:
            (init_value, (lbound, ubound), cond) = value

        lbounds[name] = lbound
        ubounds[name] = ubound
        conditions[name] = cond

        init_values[name] = init_value

    cfg.set_value('init_values', init_values)
    cfg.set_value('_lbounds', lbounds)
    cfg.set_value('_ubounds', ubounds)
    cfg.set_value('_conditions', conditions)

    # Process transformation of optimized parameters
    if cfg.get_value('transform_params'):
        from ..shared.saturation_curve import TRANSFORM_MAX_VALUE as max_value

        transform = {'ks': lambda ks: max(np.log(ks), -max_value)}
        untransform = {'ks': lambda ks_transf: min(np.exp(ks_transf), max_value)}

        SC = cfg.get_value('SC')
        if SC:
            SC.add_transformations_fns(transform, untransform,
                                       lbounds, ubounds)
    else:
        transform = untransform = None

    cfg.set_value('_transform', transform)
    cfg.set_value('_untransform', untransform)
