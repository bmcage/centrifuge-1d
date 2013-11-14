from __future__ import division
import numpy as np
from collections import OrderedDict

PARENTAL_MODULES = []

CONFIG_OPTIONS = ['inv_init_params', ('optimfn', 'leastsq'),
                  ('transform_params', True), ('untransformed_cov', False),
                  (lambda cfg: cfg.get_value('optimfn') == 'leastsq',
                    ['epsfcn', 'factor',
                     ('xtol', 1.49012e-8), ('ftol', 1.49012e-8)]),
                  (lambda cfg: (cfg.get_value('optimfn')
                                in ['fmin', 'fmin_powell']),
                    [('xtol', 1e-4), ('ftol', 1e-4), ('max_fev', None),
                     ('max_inv_iter', None), ('disp_inv_conv', True)]),
                  (lambda cfg:
                       cfg.get_value('optimfn') in ['fmin_cg', 'fmin_bfgs'],
                    [('gtol', 1e-5), ('max_inv_iter', None),
                     ('disp_inv_conv', True)])]

INTERNAL_OPTIONS = ['_transform', '_untransform', 'init_values',
                    '_lbounds', '_ubounds']

EXCLUDE_FROM_MODEL = ['inv_init_params']

def provide_inv_init_params(cfg):
    inv_init_params = cfg.get_value('inv_init_params')
    if inv_init_params and (type(inv_init_params) is dict):
        value = list(inv_init_params.keys())
    else:
        value = []

    return value

PROVIDE_OPTIONS = [provide_inv_init_params]

OPTIONS_ITERABLE_LISTS = []

def check_cfg(cfg):
    if not type(cfg.get_value('inv_init_params')) is dict:
        print("Option 'inv_init_params' has to be of type dict.")
        return False

    return True

def prior_adjust_cfg(cfg):
    pass

def adjust_cfg(cfg):
    # Process 'inv_init_params'
    init_values = OrderedDict()
    lbounds     = {}
    ubounds     = {}

    for (name, value) in cfg.get_value('inv_init_params').items():
        if value is None: continue

        (init_value, (lbound, ubound)) = value
        lbounds[name] = lbound
        ubounds[name] = ubound

        init_values[name] = init_value

    cfg.set_value('init_values', init_values)
    cfg.set_value('_lbounds', lbounds)
    cfg.set_value('_ubounds', ubounds)

    # Process transformation of optimized parameters
    if cfg.get_value('transform_params'):
        max_value = 1e150

        transform = {'ks': lambda ks: max(np.log(ks), -max_value)}
        untransform = {'ks': lambda ks_transf: min(np.exp(ks_transf), max_value)}
    else:
        transform = untransform = None

    cfg.set_value('_transform', transform)
    cfg.set_value('_untransform', untransform)
