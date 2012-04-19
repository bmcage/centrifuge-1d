PARENTAL_MODULES = []

CONFIG_OPTIONS = {
        'mandatory' : ['exp_type', 'omega', 'l0', 'l1', 'porosity',
                       'inv_init_params', 'wl_out1', 're',
                       'draw_graphs'],
        'defaults'  : {'rho': 1.0, 'g': 981.},
        'dependent' : {'theta_r':
                       (lambda cfg: len(cfg.get_value('inv_init_params')) == 2,
                        ['theta_r'])},
        'optional'  : [],
        'additional': []
        }

EXCLUDE_FROM_MODEL = []

NONITERABLE_LIST_OPTIONS = \
  (CONFIG_OPTIONS['mandatory'] + CONFIG_OPTIONS['optional'])

IGNORE_OPTIONS = []

def check_cfg(cfg):
    return True

def adjust_cfg(cfg):
    pass
