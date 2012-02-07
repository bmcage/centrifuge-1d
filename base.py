
def base_cfg():
    base = {
            'general': {'g': 981., 'debugging': False},
    'starting-filter': {'fl1': 0.0, 'ks1': -1.0 },
               'soil': {'ks': -1.0, 'l0':  -1.0, 'porosity': -1.0},
      'ending-filter': {'fl2': 0.0, 'ks2': -1.0 },
              'fluid': {'viscosity': 1.0, 'density': 1.0,
                        'wl0': 0.0, 'wl0_out': 0.0},
         'centrifuge': {'r0': -1.0, 'd0_out': -1.0,
                        'include_acceleration': True,
                        'deceleration_duration': 0.0},
         'experiment': {'exp_type': '',
                        't_start': 0.0, 't_end': 2000.0, 't_step': 200.0,
                        'omega_start': 0.0, 'omega': 35.0, 'omega_gamma': 0.5,
                        'omega_end': 0.0,
                        'data_type': 0}
    }
    return base
# unsaturated_cfg = {
#             'soil': {'n': 2.81, 'gamma': 0.0189, 'v': 1.0},
#            'fluid': {'s1_0': 0.1, 's2_0': 0.2, 'pc0': 1.0e5 },
#   'discretization': {'inner_points': 80, 'first_point_offset': 80.0,
#                      'dtype': 1, 'percent_in_saturation': 40.0,
#                      'approximation_type': 5, 'mb_epsilon': 1e-5}

def adjust_base_cfg(cfg, adjust_all = False, adjust_omega = False,
                         adjust_time = False):
    if adjust_all or adjust_omega:
        model.omega_start = model.omega_start * np.pi/ 30. # (2pi)*omega/60
        model.omega       = model.omega * np.pi/ 30.
    # TODO: adjust omega_end?

    if adjust_all or adjust_time:
        # we assure that also the last value of tspan is present (and not cut)
        model.register_key('tspan',
                           np.arange(model.t_start,
                                     model.t_end + model.deceleration_duration
                                                + model.t_step / 10.,
                                     model.t_step))


class ModelParameters:
    """
    Parameters of the centrifuge
    """
    def __init__(self, cfg = None):
        if cfg:
            self.register_keys(cfg)

    def register_key(self, key, value):
        if hasattr(self, key):
            raise Exception("Atrribute '%s' already exists !" % key)
        else:
            setattr(self, key, value)

    def register_keys(self, cfg):
        for keys_values in cfg.values():
            for (key, value) in keys_values.items():
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



