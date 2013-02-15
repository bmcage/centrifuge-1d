from __future__ import division

PARENTAL_MODULES = []

CONFIG_OPTIONS = ['inv_init_params', ('optimfn', 'leastsq'),
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

INTERNAL_OPTIONS = []

EXCLUDE_FROM_MODEL = []

PROVIDE_OPTIONS = [lambda cfg: list(cfg.get_value('inv_init_params').keys())]

OPTIONS_ITERABLE_LISTS = []

def prior_adjust_cfg(cfg):
    if not cfg.get_value('n'): cfg.set_value('n', 1.0)

def adjust_cfg(cfg):
    pass
