PARENTAL_MODULES = ['direct_draining_saturated']

CONFIG_OPTIONS = {
        'mandatory' : ['inv_init_params', 'gc1', 'wl_out1'],
        'defaults'  : {},
        'dependent' : {'inv_params':
                       (lambda cfg: len(cfg.get_value('inv_init_params')) == 2,
                        ['ks'])},
        'optional'  : [],
        'additional': []
        }

EXCLUDE_FROM_MODEL = []

IGNORE_OPTIONS = ['ks', 'n', 'gamma']

NONITERABLE_LIST_OPTIONS = ['inv_init_params']

def check_cfg(cfg):
    inv_params_len = len(cfg.get_value('inv_init_params'))
    if not ((inv_params_len == 3) or (inv_params_len == 2)):
        print('Option ''inv_params'' should be of length 3 or 2.')
        return False
    return True

def adjust_cfg(cfg):
    pass
