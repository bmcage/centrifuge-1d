PARENTAL_MODULES = []

CONFIG_OPTIONS = {
        'mandatory': ['exp_type',
                      'ks',
                      'r0', 'l0', 'duration',
                      'wt_out',
                      ],
        'defaults':  {'g': 981,
                      'omega_start': 0.0, 'omega_end': 0.0,
                      'ks1': -1.0, 'fl1': 0.0, 'ks2': -1.0, 'fl2': 0.0,
                      'fh_duration': 0.0, 'r0_fall': 1e5,
                      'density': 1.0, 'viscosity': 1.0,
                      'draw_graphs': False,
                      'include_acceleration': False,
                      'save_figures': False, 'separate_figures': False,
                      'atol': 1e-3, 'rtol':1e-6,
                      'max_steps': 800, 'max_step_size': 800.
                     },
        'dependent': {'centrifugation':
                        (lambda cfg: not (cfg.get_value('duration') == 0.0),
                         ['omega']),
                      'acceleration':
                        (lambda cfg: cfg.get_value('include_acceleration'),
                         ['omega_start', 'omega_end', 'deceleration_duration'])
                     },
        'optional':  ['l1', 'measured_gc', 'measured_w', 'gc1', 'wl_out',
                      'descr','re'],
        'additional': ['omega_fall']
        }

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
    return True

def adjust_cfg(cfg):
    """
      This method is called after the configuration is read from a file
      (and was validated). Allows to process configuration data supplied
       by configuration file(s), e.g. allocate the discretized interval
       based on the discretization type and number of inner points.
    """
    pass
