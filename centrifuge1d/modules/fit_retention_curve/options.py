from __future__ import division


from ..shared.saturation_curve import (create_SC, SC_vG)

PARENTAL_MODULES = ['base_inverse']

# These options are added by the constants.ini just as
# global defaults.ini - but as we don't need them, we
# simply ignore them. The constants.ini should probably
# depend only on the centrifuge module 'base'...
# Exception: density, g
_ignore_list =  ['tube_no', 'tube_diam', 'viscosity',
                 'ks2', 'ks1', 'fl2', 'fl1', 'fp1', 'fp2',
                 'mo_gc_calibration_curve']

def ignore_options(cfg):
    return [(name, None) for name in _ignore_list]

CONFIG_OPTIONS = ['exp_type', 'p', 'theta', 'g', 'density',
                  ('rho', 1.0),
                  'theta_s', 'theta_r',
                  ('sample_id', None), ('wl_out1', None),
                  ('measurements_filter', None),
                  ('params_ref', None),
                  ('verbosity', 1),
                  ('transform_params', False),
                  # output options
                  ('save_data', True),
                  # options generated by mkini
                  ('measurements_length', -1),
                  # ignore these defaults:
                  ignore_options
                  ]

BLACKLIST_OPTIONS = ['duration', 'deceleration_duration', 'fh_duration']

OPTIONS_ITERABLE_LISTS = []

EXCLUDE_FROM_MODEL = ['measurements_length', 'measurements_filter'] + _ignore_list

PROVIDE_OPTIONS = []

INTERNAL_OPTIONS = ['h', 'separate_figures']

def check_cfg(cfg):
    return True

def adjust_cfg(cfg):
    # Determine saturation curve model used
    cfg.set_value('sc_type', SC_vG)
    cfg.set_value('SC', create_SC(cfg))

    cfg.set_value('separate_figures', True)

    # Brute hack as default options cannot be overriden in child modules
    cfg.set_value('transform_params', False)
    cfg.set_value('_transform', None)
    cfg.set_value('_untransform', None)
