PARENTAL_MODULES = []

CONFIG_OPTIONS = ['exp_type',
                  'ks',
                  'r0', 'l0', ('l1', None), 'duration',
                  'fh_duration', 'r0_fall',
                  'wt_out',
                  'include_acceleration',
                  # solver options
                  'atol', 'rtol',
                  'first_step_size',
                  'max_steps', 'max_step_size',
                  # defaults
                  ('g', 981),
                  ('omega_start', 0.0), ('omega_end', 0.0),
                  ('ks1', -1.0), ('fl1', 0.0), ('ks2', -1.0), ('fl2', 0.0),
                  ('density', 1.0), ('viscosity', 1.0),
                  # output
                  ('draw_graphs', True),
                  ('save_figures', True), ('separate_figures', True),
                  ('save_as_text', True),
                  ('verbosity', 1),
                  # solver options
                  ('always_restart_solver', False),

                  # dependent options
                  (lambda cfg: not (cfg.get_value('duration') == 0.0),
                    ['omega']),
                  (lambda cfg: cfg.get_value('include_acceleration'),
                    ['omega_start', 'omega_end', 'deceleration_duration']),
                  ('gc1', None), ('rm1', None ),
                  ('wl0', None), ('wl1', None), ('wl_out', None),
                  ('descr', None), ('re', None)
                 ]

INTERNAL_OPTIONS = ['omega_fall', 'm']

EXCLUDE_FROM_MODEL = []

PROVIDE_OPTIONS = []

NONITERABLE_LIST_OPTIONS = []

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

    return True

def adjust_cfg(cfg):
    """
      This method is called after the configuration is read from a file
      (and was validated). Allows to process configuration data supplied
       by configuration file(s), e.g. allocate the discretized interval
       based on the discretization type and number of inner points.
    """
    from modules.shared.functions import rpm2radps, calc_omega_fall

    # Handle depending variables
    for key in ['omega', 'omega_start', 'omega_end']:
        value = cfg.get_value(key)
        if type(value) == list:
            cfg.set_value(key, [rpm2radps(omega) for omega in value])
        else:
            cfg.set_value(key, rpm2radps(value))

    cfg.set_value('omega_fall',
                  calc_omega_fall(cfg.get_value('r0_fall'), cfg.get_value('g')))


def prior_adjust_cfg(cfg):
    """
      This function is called prior to the adjust_cfg() function. It is intended
      for pre-data initialization. Mainly, if a descendent module provides an
      option, this option may not be present in the configuration, but the
      parental module may need to use the value in adjust_cfg. So here the
      needed value can be specified.
    """
    return True
