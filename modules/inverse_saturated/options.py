from __future__ import division

PARENTAL_MODULES = ['direct-saturated-heights']

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

PROVIDE_OPTIONS = ['ks']

INTERNAL_OPTIONS = ['calc_wl_out', 'calc_wl_in', 'z_size']

def adjust_cfg(cfg):
    cfg.set_value('calc_wl_in', bool(cfg.get_value('wl1')))
    cfg.set_value('calc_wl_out', bool(cfg.get_value('wl_out')))
    cfg.set_value('z_size', 2)
