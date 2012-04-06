PARENTAL_MODULES = ['base']

CONFIG_OPTIONS = {
        'mandatory' : ['wl0', 'wl1'],
        'defaults'  : {},
        'dependent' : {},
        'optional'  : [],
        'additional': ['mass_in_idx', 'mass_out_idx']
        }

EXCLUDE_FROM_MODEL = []

IGNORE_OPTIONS = []

def adjust_cfg(cfg):
    cfg.set_value('mass_in_idx',  0)
    cfg.set_value('mass_out_idx', 1)
