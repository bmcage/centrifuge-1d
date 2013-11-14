from __future__ import division

PARENTAL_MODULES = ['direct_draining_saturated', 'base_inverse']

CONFIG_OPTIONS = [(lambda cfg: cfg.get_value('optimfn') == 'raster',
                    ['raster_grid_size'])
                 ]

INTERNAL_OPTIONS = []

EXCLUDE_FROM_MODEL = []

PROVIDE_OPTIONS = []

OPTIONS_ITERABLE_LISTS = []

def check_cfg(cfg):
    from ..shared.measurements import MEASUREMENTS_NAMES
    import numpy as np

    measurements_present = False
    for name in MEASUREMENTS_NAMES.values():
        meas = cfg.get_value(name)

        if meas is None: continue

        measurements_present = True

        w_name = name + '_weights'
        weights = cfg.get_value(w_name)

        if not weights: continue

        if np.isscalar(weights):
            pass
        elif ((type(weights) == list) and
              (not (len(weights) == len(meas)))):
            print("Weights '{}'have to be a scalar or a list of the same "
                  "length as measurement.".format(w_name))
            return False
        else:
            print('Unknow weights: ' + w_name,
                  'Weights have to be a scalar or a list of the same '
                  'length as measurement.\nValue: ', weights)
            return False

    if not measurements_present:
        print('Cannot run the inverse problem as none of measured data was '
              'specified: ', MEASUREMENTS_NAMES.values())
        return False

    return True

def adjust_cfg(cfg):
    pass
