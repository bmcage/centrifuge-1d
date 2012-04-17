PARENTAL_MODULES = ['base']

def should_theta_s_be_present(cfg):
    inv_init_params_len = len(cfg.get_value('inv_init_params'))
    result = ((inv_init_params_len == 2)
              or ((inv_init_params_len == 3) and cfg.get_value('thera_r')))
    return result

def should_theta_r_be_present(cfg):
    inv_init_params_len = len(cfg.get_value('inv_init_params'))
    result = ((inv_init_params_len == 2)
              or ((inv_init_params_len == 3) and cfg.get_value('thera_r')))
    return result

CONFIG_OPTIONS = {
        'mandatory' : ['p', 'theta', 'inv_init_params'],
        'defaults'  : {'rho': 1.0},
        'dependent' : {'theta_s': (should_theta_s_be_present, ['theta_s']),
                       'theta_r': (should_theta_r_be_present, ['theta_r'])},
        'optional'  : ['sample_id'],
        'additional': []
        }

EXCLUDE_FROM_MODEL = []

NONITERABLE_LIST_OPTIONS = ['inv_init_params', 'p', 'theta']

IGNORE_OPTIONS = ['duration', 'ks', 'r0', 'l0', 'rtol', 'atol', 'wt_out', 'omega']

def check_cfg(cfg):
    theta_s_p = not cfg.get_value('thera_s') is None
    theta_r_p = not cfg.get_value('thera_r') is None

    if ((theta_s_p and theta_r_p and inv_init_params_len == 2)
        or (theta_s_p and (not theta_r_p) and inv_init_params_len == 3)
        or ((not theta_s_p) and theta_r_p and inv_init_params_len == 3)
        or ((not theta_s_p) and (not theta_r_p) and inv_init_params_len == 4)):

        # check correctness of input:
        # if 2 initial guesses given, optimize only n, gamma
        # if 3 given, we optimize also either theta_s or theta_r
        # if 4 we optimize also theta_s and theta_r
        # theta_s and theta_r;
        pass
    else:
        th_s_str = ''
        th_r_str = ''
        if theta_s_p: th_s_str = ' = ' + str(flattened_cfg['theta_s'])
        if theta_r_p: th_s_str = ' = ' + str(flattened_cfg['theta_r'])

        print('Inconsistent initial guesses inv_init_params = %s'
              % flattened_cfg['inv_init_params'])
        print("with 'theta_s'%s and 'theta_r'%s" % (th_s_str, th_r_str))
        return False

    return True

def adjust_cfg(cfg):
    pass
