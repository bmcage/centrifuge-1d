
def base_cfg():
    base = {
            'general': {'g': 981., 'debugging': False},
    'starting-filter': {'d1': 0., 'ks1': -1.0 },
               'soil': {'n': 2.81, 'gamma': 0.0189, 'ks': 2.4e-5,
                        'l': 10.0, 'porosity': 0.4, 'v': 1.0},
      'ending-filter': {'d2': 0., 'ks2': -1. },
              'fluid': {'viscosity': 1.0, 'density': 1.0,
                        's1_0': 0.1, 's2_0': 0.2, 'pc0': 1.0e5 },
         'centrifuge': {'r0': 30.0, 'l0_in': 2.0, 'l0_out': 4.0,
                        'd': 4.0, 'deceleration_duration': 0.0,
                        'include_acceleration': True},
         'experiment': {'exp_type': '',
                        't_start': 0.0, 't_end': 2000.0, 't_step': 200.0,
                        'omega_start': 0.0, 'omega': 35.0, 'omega_gamma': 0.5,
                        'omega_end': 0.0,
                        'inverse_data_filename': '', 'data_type': 0},
     'discretization': {'inner_points': 80, 'first_point_offset': 80.0, 'dtype': 1,
                        'percent_in_saturation': 40.0,
                        'approximation_type': 5, 'mb_epsilon': 1e-5}
    }
    return base


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
