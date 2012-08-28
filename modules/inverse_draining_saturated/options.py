PARENTAL_MODULES = ['direct_draining_saturated']

CONFIG_OPTIONS = ['inv_init_params',
                  # solver parameters
                  ('optimfn', 'leastsq'),
                  (lambda cfg: cfg.get_value('optimfn') == 'leastsq',
                    ['epsfcn', 'factor',
                     ('xtol', 1.49012e-8), ('ftol', 1.49012e-8)]),
                  (lambda cfg: cfg.get_value('optimfn') == 'raster',
                    ['raster_grid_size']),
                  (lambda cfg: (cfg.get_value('optimfn')
                                in ['fmin, fmin_powell']),
                    [('xtol', 1e-4), ('ftol', 1e-4)]),
                  # experiment
                  'dynamic_h_init',
                  (lambda cfg: cfg.get_value('dynamic_h_init'),
                   ['h_init_max', ('c_gammah', 1e-3)])
                 ]

INTERNAL_OPTIONS = ['calc_wl_out', 'calc_wl_in']

#EXCLUDE_FROM_MODEL = ['inv_ubounds', 'inv_lbounds']

PROVIDE_OPTIONS = [lambda cfg: list(cfg.get_value('inv_init_params').keys()),
                   'calc_gc', 'calc_rm', 'calc_wm',
                   (lambda cfg: cfg.get_value('dynamic_h_init'), ['h_init'])]

NONITERABLE_LIST_OPTIONS = ['inv_init_params', 'inv_ubounds', 'inv_lbounds',
                            'raster_grid_size']

def check_cfg(cfg):
    if ((cfg.get_value('wl_out') is None) and (cfg.get_value('gc1') is None)
        and (cfg.get_value('rm1') is None) and (cfg.get_value('wl1') is None)):
        print('No measured data (wl1, wl_out, gc1, rm1) is specified. Cannot '
              'run the inverse problem')
        return False

    return True

def prior_adjust_cfg(cfg):
    cfg.set_value('n', 1.0)

def adjust_cfg(cfg):
    from numpy import inf, asarray, power, trunc, log10

    cfg.set_value('calc_gc', bool(cfg.get_value('gc1')))
    cfg.set_value('calc_rm', bool(cfg.get_value('rm1')))
    cfg.set_value('calc_wm', False)
    cfg.set_value('calc_wl_in', bool(cfg.get_value('wl1')))
    cfg.set_value('calc_wl_out', bool(cfg.get_value('wl_out')))
