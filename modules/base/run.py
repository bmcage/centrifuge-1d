from os.path import exists
import numpy as np

PARAMETERS =  {
            'general': {'g': 981., 'debugging': False},
    'starting-filter': {'fl1': 0.0, 'ks1': -1.0 },
    'soil': {#'ks': -1.0, 'l0':  -1.0, 'porosity': -1.0
             },
      'ending-filter': {'fl2': 0.0, 'ks2': -1.0 },
              'fluid': {'viscosity': 1.0, 'density': 1.0,
    #                        'wl0': 0.0, 'wl0_out': 0.0
                        },
    'centrifuge': {#'r0': -1.0, 'r0_fall': 1.0e5, 'd0_out': -1.0,
                        'include_acceleration': True,
                        'deceleration_duration': 0.0},
         'experiment': {'exp_type': '',
    #                        't_start': 0.0, 't_end': 2000.0, 't_step': 200.0,
    #                    'omega_start': -1.0, 'omega': [35.0], 'omega_gamma': 0.5,
    #                    'omega_end': -1.0, 
                        'duration': -1.0,
                        'data_type': 0}
    }

def base_cfg():
    base = {
         'general': {'g': 981. }}
    return base

# unsaturated_cfg = {
#             'soil': {'n': 2.81, 'gamma': 0.0189, 'v': 1.0},
#            'fluid': {'s1_0': 0.1, 's2_0': 0.2, 'pc0': 1.0e5 },
#   'discretization': {'inner_points': 80, 'first_point_offset': 80.0,
#                      'dtype': 1, 'percent_in_saturation': 40.0,
#                      'approximation_type': 5, 'mb_epsilon': 1e-5}

def check_attributes(cfg, attributes):
    for attribute in attributes:
        if not attribute in cfg:
            print('Mandatory argument is not present: ', attribute)
            exit(1)

def adjust_cfg(cfg):

    mandatory_attributes = ['omega', 'duration']

    check_attributes(cfg, mandatory_attributes)
    if cfg['include_acceleration']:
        check_attributes(cfg, ['omega_start', 'omega_end'])

    omega_fn = lambda x: x * np.pi/ 30.0 # omega -> (2pi)*omega/60//rpm->rad.s-1

    for attribute in ['omega', 'omega_start', 'omega_end']:
        if attribute in cfg:
            value = cfg[attribute]
            if np.isscalar(value):
                cfg[attribute] = omega_fn(value)
            else:
                cfg[attribute] = [omega_fn(omega) for omega in value]

    # we assure that also the last value of tspan is present (and not cut)
    if 'duration' in cfg:
        if ('t_start' in cfg) or ('t_end' in cfg) or ('t_step' in cfg):
            raise ValueError('Cannot set both: "duration" and one of "t_start",'
                             ' "t_end", t_step"')
    else:
        cfg['duration'] = (cfg['t_end'] - cfg['t_start']) / cfg['t_step']

def generate_tubes_suffixes(tubes_numbers):
    if not tubes_numbers:
        print('CFG: Error: No tubes numbers specified: "%s"' % tubes_numbers)
        exit(1)

    suffixes = ['-tube' + str(tube_no) for tube_no in tubes_numbers]
    identifiers = [', tube ' + str(tube_no) for tube_no in tubes_numbers]

    return suffixes, identifiers
