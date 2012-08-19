PARENTAL_MODULES = ['base']

CONFIG_OPTIONS = ['wl0', 'wl1'],

INTERNAL_OPTIONS = ['mass_in_idx', 'mass_out_idx']

def adjust_cfg(cfg):
    cfg.set_value('mass_in_idx',  0)
    cfg.set_value('mass_out_idx', 1)
