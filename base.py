import numpy as np

def base_cfg():
    base = {
            'general': {'g': 981., 'debugging': False},
    'starting-filter': {'fl1': 0.0, 'ks1': -1.0 },
               'soil': {'ks': -1.0, 'l0':  -1.0, 'porosity': -1.0},
      'ending-filter': {'fl2': 0.0, 'ks2': -1.0 },
              'fluid': {'viscosity': 1.0, 'density': 1.0,
                        'wl0': 0.0, 'wl0_out': 0.0},
         'centrifuge': {'r0': -1.0, 'r0_fall', 'd0_out': -1.0,
                        'include_acceleration': False,
                        'deceleration_duration': 0.0},
         'experiment': {'exp_type': '',
                        't_start': 0.0, 't_end': 2000.0, 't_step': 200.0,
                        'omega_start': -1.0, 'omega': 35.0, 'omega_gamma': 0.5,
                        'omega_end': -1.0, 'duration': -1.0,
                        'data_type': 0}
    }
    return base
# unsaturated_cfg = {
#             'soil': {'n': 2.81, 'gamma': 0.0189, 'v': 1.0},
#            'fluid': {'s1_0': 0.1, 's2_0': 0.2, 'pc0': 1.0e5 },
#   'discretization': {'inner_points': 80, 'first_point_offset': 80.0,
#                      'dtype': 1, 'percent_in_saturation': 40.0,
#                      'approximation_type': 5, 'mb_epsilon': 1e-5}

def adjust_cfg(cfg):
    # print(cfg)
    omega_fn = lambda x: x * np.pi/ 30.0 # omega -> (2pi)*omega/60//rpm->rad.s-1

    cfg['omega'] = map(omega_fn, cfg['omega'])
    if cfg['include_acceleration']:
        if 'omega_start' in cfg and not (cfg['omega_start'] < 0.):
            cfg['omega_start'] =  map(omega_fn, cfg['omega_start']) # (2pi)*omega/60
        else:
            cfg['omega_start'] = np.asarray([0.0], dtype=float)

        if 'omega_end' in cfg and not (cfg['omega_end'] < 0.):
            cfg['omega_end'] = map(omega_fn, cfg['omega_end']) # (2pi)*omega/60
        else:
            cfg['omega_end'] = np.asarray([0.0], dtype=float)
    else:
        cfg['omega_start'] = cfg['omega']
        cfg['omega_end']   = cfg['omega']

    # we assure that also the last value of tspan is present (and not cut)
    if 'duration' in cfg:
        if ('t_start' in cfg) or ('t_end' in cfg) or ('t_step' in cfg):
            raise ValueError('Cannot set both: "duration" and one of "t_start",'
                             ' "t_end", t_step"')
    else:
        cfg['duration'] = (cfg['t_end'] - cfg['t_start']) / cfg['t_step']

class ModelParameters:
    """
    Parameters of the centrifuge
    """
    def __init__(self, cfg = None):
        if cfg:
            self.register_keys(cfg)
        self.register_key('tspan', np.asarray([], dtype=float))
        self.register_key('omega_fall', np.asarray([], dtype=float))

    def register_key(self, key, value):
        if hasattr(self, key):
            raise Exception("Atrribute '%s' already exists !" % key)
        else:
            if type(value) == str:
                setattr(self, key, value)
            else:
                setattr(self, key, np.asarray(value))

    def register_keys(self, flattened_cfg):
        for (key, value) in flattened_cfg.items():
            self.register_key(key.lower(), value)

    def set(self, key, value):
        key_lower = key.lower()

        if not hasattr(self, key_lower):
            raise Exception('Atrribute ''%s'' does not exist !' % key)

        value_type = type(value)
        key_type   = type(getattr(self, key_lower))

        if value_type == key_type:
            setattr(self, key_lower, value)
        elif value_type == int and key_type == float:
            #convert int automatically to float
            setattr(self, key, float(value))
        elif value_type == list:
            for item in value:
                if type(item) == int and key_type == float:
                    pass
                elif not ((type(item) == key_type) or
                          (type(item) == int and key_type == float)):
                    raise ValueError("ModelParameters: key '%s' has wrong type."
                            " Expected type '%s' and got type '%s' of %s"
                            % (key, key_type, value_type, value))
                if value and type(value[0] == int) and (key_type == float):
                    value = [float(item) for item in value]

                setattr(self, key, value)
        else:
            raise ValueError("ModelParameters: key '%s' has wrong type. Expected"
                             " type '%s' and got type '%s' of %s"
                             % (key, key_type, value_type, value))



