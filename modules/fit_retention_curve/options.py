PARENTAL_MODULES = ['base']

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

CONFIG_OPTIONS = {
        'mandatory' : ['p', 'theta', 'inv_init_params'],
        'defaults'  : {'rho': 1.0},
        'dependent' : {'theta_s': (should_theta_s_be_present, ['theta_s']),
                       'theta_r': (should_theta_r_be_present, ['theta_r'])},
        'optional'  : ['sample_id', 'wl_out1'],
        'additional': []
        }

EXCLUDE_FROM_MODEL = []

NONITERABLE_LIST_OPTIONS = ['inv_init_params', 'p', 'theta']

PROVIDE_OPTIONS = ['duration', 'ks', 'r0', 'l0', 'wt_out',
                  'omega']

BLACKLIST_OPTIONS = ['rtol', 'atol', 'include_acceleration', 'max_step_size',
                     'r0_fall', 'fh_duration', 'max_steps']

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

def adjust_cfg(cfg):
    pass
