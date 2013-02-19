from __future__ import division

PARENTAL_MODULES = ['base']

CONFIG_OPTIONS = [('porosity', None)]

INTERNAL_OPTIONS = ['mass_in_idx', 'mass_out_idx', 'z_size']

OPTIONS_ITERABLE_LISTS = ['l0', 'wl0']

def adjust_cfg(cfg):
    cfg.set_value('mass_in_idx',  0)
    cfg.set_value('mass_out_idx', 1)
    cfg.set_value('z_size', 2)

    cfg.set_value('calc_gf_mt', False)

def check_cfg(cfg):
    if not (cfg.get_value('wl0') or cfg.get_value('ww0')):
        print("One of 'wl0' or 'ww0' parameters must be specified.")
        return False

    return True
