import numpy as np
from modules.shared.utils import test_related_options
from modules.shared.measurements import MEASUREMENTS_NAMES

PARENTAL_MODULES = []

def gen_measurements(cfg):
    names   = [(name, None) for name in MEASUREMENTS_NAMES.values() \
               if not name is None]
    xvalues = [(name + '_xvalues', None) for name in MEASUREMENTS_NAMES.values()  \
               if not name is None]
    weights = [(name + '_weights', 1.0) for name in MEASUREMENTS_NAMES.values() \
               if not name is None]

    return names + xvalues + weights

CONFIG_OPTIONS = ['exp_type',
                  'ks',
                  ('r0', None), ('re', None), ('wl0', None), 'l0', 'duration',
                  ('ww0', None),
                  'fh_duration', 'r0_fall',
                  ('acceleration_duration', 21.),
                  'include_acceleration',
                  # solver options
                  'atol', 'rtol',
                  'first_step_size',
                  'max_steps', 'max_step_size',
                  # defaults
                  ('g', 981),
                  ('tube_no', None),
                  ('tube_diam', None),
                  'fl1',
                  (lambda cfg: cfg.get_value('fl1') > 0.0,
                      ['ks1'], [('ks1', -1.0)]),
                  'fl2',
                  (lambda cfg: cfg.get_value('fl2') > 0.0,
                      ['ks2'], [('ks2', -1.0)]),
                  ('density', 1.0), ('viscosity', 1.0),
                  # output
                  ('show_figures', True),
                  ('save_figures', False), ('separate_figures', False),
                  ('save_data', True),
                  ('verbosity', 1),
                  # solver options
                  ('always_restart_solver', False),

                  # dependent options
                  (lambda cfg: not (cfg.get_value('duration') == 0.0),
                    ['omega']),
                  (lambda cfg: cfg.get_value('include_acceleration'),
                    ['deceleration_duration'],
                    [('deceleration_duration', 0.0)]),
                  # measurements and referencing parameters
                  ('smoothing', {}),
                  gen_measurements,
                  ('ww1', None), ('gf_mo_tara', None), ('gf_mt_tara', None),
                  ('measurements_scale_coefs', None),
                  ('cf_weights', 1.0),
                  ('measurements_keep', None), ('measurements_remove', None),

                  ('descr', None), ('re', None),
                  (lambda cfg: not (cfg.get_value('gf_mo') is None),
                    [('gf_mo_calibration_curve', None)]),
                  (lambda cfg: not (cfg.get_value('gf_mt') is None),
                    [('gf_mt_calibration_curve', None)]),
                  ('mo_gc_calibration_curve', None),

                  # options generated by mkini
                  ('measurements_length', -1)
                 ]

INTERNAL_OPTIONS = ['omega2g_fns', 'find_omega2g', 't0', 'omega_start',
                    'calc_gf_mt', 'calc_gf_mo', 'calc_gcf_mo',
                    'tube_crosssectional_area', 'z_size']

EXCLUDE_FROM_MODEL = ['measurements_length', 'omega2g_fns', 'r0',
                      'gf_mt_calibration_curve', 'gf_mo_calibration_curve',
                      'tube_diam', 'ww0', 'ww1', 'gf_mo_tara', 'gf_mt_tara']

PROVIDE_OPTIONS = []

OPTIONS_ITERABLE_LISTS = ['r0', 're', 'l0', 'omega', 'porosity',
                           'duration', 'fh_duration', 'deceleration_duration']

def check_cfg(cfg):
    """
      Perform additional test on 'cfg' to validate it. Test for expected value,
      type etc. of supplied values should be done here. Checking for the
      presence of mandatory and dependent options  is done by default.
      Additional information informing user about failed test(s) should be
      supplied. The return value of this function is type boolean.
    """
    if not cfg.get_value('include_acceleration'):
        value = cfg.get_value('deceleration_duration')
        is_list = (type(value) in [list, tuple])
        if (is_list and any(value)) or (not is_list and value):
            print("Option 'deceleration_duration' can't have a positive value "
                  "if 'include_acceleration' is False.")
            return False

    for options_names in (('r0', 're'), ):
        if not test_related_options(cfg, options_names, mode='exact1'):
            return False

    for options_names in (('wl0', 'ww0'), ('wl1', 'ww1')):
        if not test_related_options(cfg, options_names, mode='atmost1'):
            return False

    if cfg.get_value('tube_diam') is None:
        print('Tube diameter is not specified.')
        return False

    for F_name in ('gf_mo', 'gf_mt'):
        F = cfg.get_value(F_name)

        if (F is None): continue

        F_tara = cfg.get_value(F_name + '_tara')
        if not F_tara is None:
            if ((not type(F_tara) in (list, tuple))
                or (len(F_tara) != 2)):
                print('The tara value of ' + F_name + ' (' + F_name + '_tara) '
                      'must be a list/tuple of length 2 of type '
                      '[omega, weight]. Aborting.')
                return False

        F_calibration_curve = cfg.get_value(F_name + '_calibration_curve')

        if type(F_calibration_curve) in (list, tuple):
            if not (len(F_calibration_curve) == len(F)):
                print("Force calibration curve '" + F_name + "_calibration_curve"
                      "' supplied as array has to be of the same length "
                    "as the measured force '"+ F + "'")
                return False
        elif type(F_calibration_curve) in (float, int, dict):
            pass
        elif F_calibration_curve is None:
            print('INFO: Calibration curve: ' +  F_name + '_calibration_curve'
                  + ' was not specified, difference between two subsequent '
                  'forces will be used instead.')
            pass
        else:
            print('Unsuppported type for calibration_curve of ' + F_name + '.'
                  '\nOnly float/int/array of floats/dict or ints is allowed.')
            return False

        if F_name in ('gf_mo'):
            MO_GC = cfg.get_value('mo_gc_calibration_curve')

            if MO_GC is None:
                print('No calibration curve for expelled water was '
                      'specified. Cannot continue. Exiting...')
                return False

            if ((not type(MO_GC) in (list, tuple))
                or (not type(MO_GC[0]) in (list, tuple))):
                print("Calibration curve must by of type 'list' or 'tuple': ",
                      MO_GC)
                return False

            if not len(MO_GC) > 1:
                print('Calibration curve must by of length at least 2: ',
                      MO_GC)
                return False

        if not cfg.get_value('mo_gc_calibration_curve'):
            print("WARNING: 'MO_GC_calibration_curve' is not set. This may "
                  "lead to problems. You have been warned:)")

    return True

def adjust_cfg(cfg):
    """
      This method is called after the configuration is read from a file
      (and was validated). Allows to process configuration data supplied
       by configuration file(s), e.g. allocate the discretized interval
       based on the discretization type and number of inner points.
    """
    from modules.shared.functions import (find_omega2g, find_omega2g_fh,
                                          find_omega2g_dec)

    # NOTE: conversion of omega from rpm to radians is done inside the
    #       read() method of class Measurements
    cfg.set_value('omega2g_fns', {'a': find_omega2g, 'g': find_omega2g_fh,
                                  'd': find_omega2g_dec})

    # if r0 was set (but not rE), we set rE (and discard r0)
    if not cfg.get_value('re'):
        r0_np  = np.asarray(cfg.get_value('r0'), dtype=float)
        l0_np  = np.asarray(cfg.get_value('l0'), dtype=float)
        fl2_np = np.asarray(cfg.get_value('fl2'), dtype=float)

        rE = r0_np + l0_np + fl2_np
        if not np.isscalar(rE):
            rE = list(rE)
        cfg.set_value('re', rE)

    cfg.set_value('calc_gf_mt', True)
    cfg.set_value('calc_gcf_mo', False)
    if not cfg.get_value('mo_gc_calibration_curve') is None:
        cfg.set_value('calc_gf_mo', True)

    tube_crosssectional_area = \
      np.pi * np.power(cfg.get_value('tube_diam'), 2) / 4.
    cfg.set_value('tube_crosssectional_area', tube_crosssectional_area)

    # if ww0 (ww1) was provided, we set wl0 (wl1)
    for name in ('ww0', 'ww1'):
        value = cfg.get_value(name)

        if value:
            new_value = (np.asarray(value, dtype=float)
                         / cfg.get_value('density') / tube_crosssectional_area)
            if np.isscalar(new_value):
                cfg.set_value('wl'+name[2], [new_value])
            else:
                cfg.set_value('wl'+name[2], list(new_value))

def prior_adjust_cfg(cfg):
    """
      This function is called prior to the adjust_cfg() function. It is intended
      for pre-data initialization. Mainly, if a descendent module provides an
      option, this option may not be present in the configuration, but the
      parental module may need to use the value in adjust_cfg. So here the
      needed value can be specified.
    """
    return True
