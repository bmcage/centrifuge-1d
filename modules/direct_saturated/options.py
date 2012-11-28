from __future__ import division

PARENTAL_MODULES = ['base']

CONFIG_OPTIONS = [('porosity', None)]

INTERNAL_OPTIONS = ['mass_in_idx', 'mass_out_idx']

OPTIONS_ITERABLE_LISTS = ['l0', 'wl0']

def adjust_cfg(cfg):
    cfg.set_value('mass_in_idx',  0)
    cfg.set_value('mass_out_idx', 1)

def check_cfg(cfg):
    if not cfg.get_value('wl0'):
        print('Parameter \'wl0\' is not specified.')
        return False

    return True
