from __future__ import division

PARENTAL_MODULES = ['direct_draining_saturated']

CONFIG_OPTIONS = ['inv_init_params',
                  # solver parameters
                  ('optimfn', 'leastsq'),
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
                     ('disp_inv_conv', True)]),
                  (lambda cfg: cfg.get_value('optimfn') == 'raster',
                    ['raster_grid_size']),
                  # measurement weights
                  ('wl1_weights', None), ('wl_out_weights', None),
                  ('gc1_weights', None), ('rm1_weights', None),
                  ('cf_weights', None)
                 ]

INTERNAL_OPTIONS = []

#EXCLUDE_FROM_MODEL = ['inv_ubounds', 'inv_lbounds']

PROVIDE_OPTIONS = [lambda cfg: list(cfg.get_value('inv_init_params').keys())]

OPTIONS_ITERABLE_LISTS = []

def check_cfg(cfg):
    from modules.shared.characteristics import MEASUREMENTS_NAMES

    measurements_present = False
    for name in MEASUREMENTS_NAMES.values():
        meas = cfg.get_value(name)

        if meas is None: continue

        measurements_present = True

        w_name = name + '_weights'
        weights = cfg.get_value(w_name)

        if not weights: continue

        if not type(weights) == list:
            print('Weights has to be a list: {}'.format(w_name))
            return False

        if not (len(weights) == len(meas)):
            print('Weights have to be of the same length as measurement: '
                  '{}'.format(name))
            return False

    if not measurements_present:
        print('Cannot run the inverse problem as none of measured data was '
              'specified: ', MEASUREMENTS_NAMES.values())
        return False

    return True

def prior_adjust_cfg(cfg):
    cfg.set_value('n', 1.0)

def adjust_cfg(cfg):
    pass
