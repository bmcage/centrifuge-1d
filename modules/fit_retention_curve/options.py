PARENTAL_MODULES = ['base']

CONFIG_OPTIONS = {
        'mandatory' : ['p', 'theta', 'inv_init_params'],
        'defaults'  : {'rho': 1.0},
        'dependent' : {},
        'optional'  : ['sample_id'],
        'additional': []
        }

EXCLUDE_FROM_MODEL = []

NONITERABLE_LIST_OPTIONS = ['inv_init_params', 'p', 'theta']

IGNORE_OPTIONS = ['duration', 'ks', 'r0', 'l0', 'rtol', 'atol', 'wt_out', 'omega']

def adjust_cfg(cfg):
    pass