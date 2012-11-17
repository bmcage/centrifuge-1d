from __future__ import division

PARENTAL_MODULES = ['base', 'direct_draining_saturated']

CONFIG_OPTIONS = ['percent_in_saturation', 30]

INTERNAL_OPTIONS = ['saturated_idx']

EXCLUDE_FROM_MODEL = ['percent_in_saturation']

def check_cfg(cfg):
    value = cfg.get_value('percent_in_saturation')

    if ((not type(value) == int) and (not type(value) == float)
        or (value > 90) or (value < 10)):

        print('cfg_check error: the value of ''percent_in_saturation'' has '
              'to be between 10 and 90')
        return False

    return True

def adjust_cfg(cfg):
    cfg.set_value('saturated_idx',
                  int(cfg.get_value('inner_points') / 100.
                      * (100. - float(cfg.get_value('percent_in_saturation')))))
