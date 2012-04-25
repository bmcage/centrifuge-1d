PARENTAL_MODULES = ['direct_draining_saturated']

CONFIG_OPTIONS = {
        'mandatory' : ['inv_init_params'],
        'defaults'  : {},
        'dependent' : {'inv_params':
                       (lambda cfg: len(cfg.get_value('inv_init_params')) == 2,
                        ['ks'])},
        'optional'  : ['wl_out1', 'gc1', 'rm1'],
        'additional': ['calc_wl_out']
        }

EXCLUDE_FROM_MODEL = []

PROVIDE_OPTIONS = ['ks', 'n', 'gamma', 'calc_gc', 'calc_rm']

NONITERABLE_LIST_OPTIONS = ['inv_init_params']

def check_cfg(cfg):
    inv_params_len = len(cfg.get_value('inv_init_params'))
    if not ((inv_params_len == 3) or (inv_params_len == 2)):
        print('Option ''inv_params'' should be of length 3 or 2.')
        return False

    if ((cfg.get_value('wl_out1') is None) and (cfg.get_value('gc1') is None)
        and (cfg.get_value('rm1') is None)):
        print('No measured data (wl_out1, gc1, rm1) is specified. Cannot run '
              'inverse problem')
        return False

    return True

def adjust_cfg(cfg):
    cfg.set_value('calc_gc', not cfg.get_value('gc1') is None)
    cfg.set_value('calc_rm', not cfg.get_value('rm1') is None)
    cfg.set_value('calc_wl_out', not cfg.get_value('wl_out1') is None)
