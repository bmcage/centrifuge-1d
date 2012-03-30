PARENTAL_MODULES = []

CONFIG_OPTIONS = {
        'mandatory': ['exp_type',
                      'r0', 'l0', 'omega', 'duration',
                      'include_acceleration', ,
                      'wt_out',
                      'atol', 'rtol'],
        'defaults': {'g': 981,
                     'omega_start': 0.0, 'omega_end': 0.0,
                     'ks1': -1.0, 'fl1': 0.0, 'ks2': -1.0, 'fl2': 0.0,
                     'fh_duration': 0.0, 'r0_fall': 1e5,
                     'density': 1.0, 'viscosity': 1.0},
        'dependent': {'acceleration': \
                      (lambda cfg: cfg.get_value('include_acceleration'),
                       ['omega_start', 'omega_end', 'deceleration_duration'])},
        'optional': ['l1']
        'additional': []
        }

CONFIG_ONLY_OPTIONS = ['exp_type']

IGNORE_OPTIONS = []
