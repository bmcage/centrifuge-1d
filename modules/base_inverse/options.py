from __future__ import division
import numpy as np

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

INTERNAL_OPTIONS = ['_transform', '_untransform']

EXCLUDE_FROM_MODEL = []

PROVIDE_OPTIONS = [lambda cfg: list(cfg.get_value('inv_init_params').keys())]

OPTIONS_ITERABLE_LISTS = []

def prior_adjust_cfg(cfg):
    pass

def adjust_cfg(cfg):
    if cfg.get_value('transform_params'):
        max_value = 1e150

        transform = {'ks': lambda ks: max(np.log(ks), -max_value),
                     'n':  lambda n: max(np.log(n - 1.0), -max_value),
                     'gamma': lambda gamma: max(np.log(-gamma), -max_value)}
        untransform = {'ks': lambda ks_transf: min(np.exp(ks_transf), max_value),
                       'n': lambda n_transf: 1+min(np.exp(n_transf), max_value),
                       'gamma': lambda gamma_transf: -min(np.exp(gamma_transf), max_value)}
    else:
        transform = untransform = None

    cfg.set_value('_transform', transform)
    cfg.set_value('_untransform', untransform)
