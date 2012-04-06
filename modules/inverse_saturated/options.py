PARENTAL_MODULES = ['direct-saturated-heights']

CONFIG_OPTIONS = {
        'mandatory' : ['inv_init_params'],
        'defaults'  : {},
        'dependent' : {},
        'optional'  : [],
        'additional': ['mass_in_idx', 'mass_out_idx']
        }

EXCLUDE_FROM_MODEL = []

IGNORE_OPTIONS = ['ks']

def adjust_cfg(cfg):
    pass
