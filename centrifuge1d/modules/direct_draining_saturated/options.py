from __future__ import division

PARENTAL_MODULES = ['base_unsaturated']

CONFIG_OPTIONS = [('dynamic_h_init', False),
                  (lambda cfg: cfg.get_value('dynamic_h_init'),
                   ['h_init_max', ('c_gammah', 1e-3)], ['h_init']),
                  'rb_type',
                  # dependent
                  (lambda cfg: cfg.get_value('rb_type') == 2,
                    ['h_last']),
                  (lambda cfg: cfg.get_value('rb_type') == 3,
                    ['dip_height', 'h_last',
                     ('s2_0', None), ('s2_atol', 1e-8)]),
                  'h_last'
                 ]

INTERNAL_OPTIONS = ['first_idx', 'last_idx', 's1_idx', 's2_idx',
                    'mass_in_idx', 'mass_out_idx', 'pq_idx',
                    'wm0', 'z_size']

EXCLUDE_FROM_MODEL = []

PROVIDE_OPTIONS = []

OPTIONS_ITERABLE_LISTS = []

def adjust_cfg(cfg):
    # Set array indexes
    inner_points = cfg.get_value('inner_points')
    cfg.set_value('sc_max_refine', 0)
    cfg.set_value('first_idx',    0)
    cfg.set_value('last_idx',     inner_points+1)
    cfg.set_value('mass_in_idx',  inner_points+2)
    cfg.set_value('s1_idx',       inner_points+3)
    cfg.set_value('s2_idx',       inner_points+4)
    cfg.set_value('mass_out_idx', inner_points+5)
    cfg.set_value('pq_idx',       inner_points+6)
    # total length of 'z' array (discretization points + s1,s2,mass_in,...)
    cfg.set_value('z_size',       inner_points+7)

    # if 'dynamic_h_init' then 'h_init' is set at runtime
