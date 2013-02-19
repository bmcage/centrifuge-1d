from __future__ import division
import modules.shared.saturation_curve as mSC

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
                  ('rho', 1.0), 'show_figures',
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

OPTIONS_ITERABLE_LISTS = []

EXCLUDE_FROM_MODEL = ['measurements_length', 'measurements_filter'] + _ignore_list

PROVIDE_OPTIONS = []

INTERNAL_OPTIONS = ['h', 'separate_figures']

def check_cfg(cfg):
    return True

def filter_measurements(cfg, measurements):
    fltr = cfg.get_value('measurements_filter')

    if not fltr: return

    for (name, value) in cfg.iterate_values():
        if (not type(value) in [list, tuple]) or (not name in measurements):
             continue

        filtered_value = []
        for (v, keep_p) in zip(value, fltr):
            if keep_p: filtered_value.append(v)

        print('name:', name, '\nvalue:', value, 'fv:', filtered_value,
              '\nfilter:', fltr)
        cfg.set_value(name, filtered_value)

def adjust_cfg(cfg):
    filter_measurements(cfg, ['theta', 'p'])
    cfg.set_value('separate_figures', True)

    # Brute hack as default options cannot be overriden in child modules
    cfg.set_value('transform_params', False)
    cfg.set_value('_transform', None)
    cfg.set_value('_untransform', None)

    SC = mSC.SC_vanGenuchten()
    cfg.set_value('SC', SC)
