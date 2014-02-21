from __future__ import division, print_function
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
                    '_lbounds', '_ubounds', '_conditions']

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

    SC_type = cfg.get_value('sc_type')
    if (SC_type == 1):
        n = init_values.get('n') or cfg.get_value('n')
        gamma = init_values.get('gamma') or cfg.get_value('gamma')
        if not n or not gamma:
            print("Option SC_type = 1 requires init 'n' and 'gamma'")
            exit(1)
    elif (SC_type == 2):
        hi = init_values.get('hi') or cfg.get_value('hi')
        ki = init_values.get('ki') or cfg.get_value('ki')
        ui = init_values.get('ui') or cfg.get_value('ui')
        if not hi or not ki or not ui:
            print("Option SC_type = 2 requires init 'hi',  'ki' and 'ui'")
            exit(1)
        if not len(hi) == len(ki) == len(ui):
            print ("Length of hi, ki and ui must be identical for SC_type = 2")
            exit(1)

    # Process transformation of optimized parameters
    if cfg.get_value('transform_params'):
        max_value = 1e150

        transform = {'ks': lambda ks: max(np.log(ks), -max_value)}
        untransform = {'ks': lambda ks_transf: min(np.exp(ks_transf), max_value)}

        SC = cfg.get_value('SC')
        if SC:
            SC.add_transformations_fns(transform, untransform, max_value)
    else:
        transform = untransform = None

    cfg.set_value('_transform', transform)
    cfg.set_value('_untransform', untransform)
