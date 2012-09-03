PARENTAL_MODULES = ['direct-saturated-heights']

CONFIG_OPTIONS = ['inv_init_params']

PROVIDE_OPTIONS = ['ks']

INTERNAL_OPTIONS = ['calc_wl_out', 'calc_wl_in']

def adjust_cfg(cfg):
    cfg.set_value('calc_wl_in', bool(cfg.get_value('wl1')))
    cfg.set_value('calc_wl_out', bool(cfg.get_value('wl_out')))
