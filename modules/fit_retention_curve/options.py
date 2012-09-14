PARENTAL_MODULES = []

def should_theta_s_be_present(cfg):
    inv_init_params_len = len(cfg.get_value('inv_init_params'))
    result = ((inv_init_params_len == 2)
              or ((inv_init_params_len == 3)
                   and (cfg.get_value('theta_r') is None)))
    return result

def should_theta_r_be_present(cfg):
    inv_init_params_len = len(cfg.get_value('inv_init_params'))
    result = ((inv_init_params_len == 2)
              or ((inv_init_params_len == 3)
                  and (cfg.get_value('theta_s') is None)))
    return result

CONFIG_OPTIONS = ['exp_type', 'g', 'p', 'theta', 'inv_init_params',
                  ('rho', 1.0), 'show_figures',
                  (should_theta_s_be_present, ['theta_s']),
                  (should_theta_r_be_present, ['theta_r']),
                  ('sample_id', None), ('wl_out1', None),
                  ('measurements_filter', None),
                  ('n_ref', None), ('gamma_ref', None),

                  # options generated by mkini
                  ('measurements_length', -1)]

NONITERABLE_LIST_OPTIONS = ['inv_init_params', 'p', 'theta',
                            'measurements_filter']

EXCLUDE_FROM_MODEL = ['measurements_length', 'measurements_filter']

INTERNAL_OPTIONS = []

def check_cfg(cfg):
    theta_s_p = not cfg.get_value('theta_s') is None
    theta_r_p = not cfg.get_value('theta_r') is None

    inv_init_params_len = len(cfg.get_value('inv_init_params'))

    # check correctness of input:
    # if 2 initial guesses (parameters) are given, optimize only n, gamma
    # if 3 are given, we optimize also either theta_s or theta_r
    # if 4 we optimize also both theta_s and theta_r
    if not ((theta_s_p and theta_r_p and inv_init_params_len == 2)
            or (theta_s_p and (not theta_r_p) and inv_init_params_len == 3)
            or ((not theta_s_p) and theta_r_p and inv_init_params_len == 3)
            or ((not theta_s_p) and (not theta_r_p)
                and inv_init_params_len == 4)):
        th_s_str = ''
        th_r_str = ''
        if theta_s_p: th_s_str = ' = ' + str(cfg.get_value('theta_s'))
        if theta_r_p: th_s_str = ' = ' + str(cfg.get_value('theta_r'))

        print('Inconsistent initial guesses inv_init_params = %s'
              % cfg.get_value('inv_init_params'))
        print("with 'theta_s'%s and 'theta_r'%s" % (th_s_str, th_r_str))
        return False

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
