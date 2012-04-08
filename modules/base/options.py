PARENTAL_MODULES = []

CONFIG_OPTIONS = {
        'mandatory': ['exp_type',
                      'ks',
                      'r0', 'l0', 'duration',
                      'wt_out',
                      'atol', 'rtol'],
        'defaults':  {'g': 981,
                      'omega_start': 0.0, 'omega_end': 0.0,
                      'ks1': -1.0, 'fl1': 0.0, 'ks2': -1.0, 'fl2': 0.0,
                      'fh_duration': 0.0, 'r0_fall': 1e5,
                      'density': 1.0, 'viscosity': 1.0,
                      'draw_graphs': False,
                      'include_acceleration': False
                     },
        'dependent': {'centrifugation':
                        (lambda cfg: not (cfg.get_value('duration') == 0.0),
                         ['omega']),
                      'acceleration':
                        (lambda cfg: cfg.get_value('include_acceleration'),
                         ['omega_start', 'omega_end', 'deceleration_duration'])
                     },
        'optional':  ['l1', 'measured_gc', 'measured_w'],
        'additional': ['omega_fall']
        }

EXCLUDE_FROM_MODEL = []

IGNORE_OPTIONS = []

def adjust_cfg(cfg):
    """
      This method is called after the configuration is read from a file
      (and was validated). Allows to process configuration data supplied
       by configuration file(s), e.g. allocate the discretized interval
       based on the discretization type and number of inner points.
    """
    duration = cfg.get_value('duration')
    if (type(duration) == list) and (not any(duration)) or (duration == 0.0):
        # We perform falling-head test only
        cfg.set_value('duration', cfg.get_value('fh_duration'))
        cfg.set_value('fh_duration', 0.0)

        cfg.set_value('include_acceleration', False)
        cfg.set_value('r0', cfg.get_value('r0_fall'))
        cfg.set_value('omega', cfg.get_value('omega_fall'))
