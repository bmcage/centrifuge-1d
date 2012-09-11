from modules.shared.functions import lagrangean_derivative_coefs
from numpy import linspace, power, asarray

PARENTAL_MODULES = ['base']

CONFIG_OPTIONS = ['inner_points', 'dtype', 'n', 'gamma',
                  'h_init', 'porosity',
                  'calc_gc', 'calc_rm', 'calc_wm',
                  'rb_type',
                  # dependent
                  (lambda cfg: cfg.get_value('rb_type') == 2,
                    ['h_last']),
                  (lambda cfg: cfg.get_value('rb_type') == 3,
                    ['dip_height', 'h_last',
                     ('s2_0', None), ('s2_atol', 1e-8)]),
                  # optional
                  ('n1', None), ('gamma1', None),
                  ('n2', None), ('gamma2', None),
                  'h_last'
                 ]

INTERNAL_OPTIONS = ['m', 'y', 'y12', 'dy', 'ldc1', 'ldc2', 'ldc3',
                    'first_idx', 'last_idx', 's1_idx', 's2_idx',
                    'mass_in_idx', 'mass_out_idx', 'pq_idx', 'z_size',
                    'wm0']

EXCLUDE_FROM_MODEL = ['dtype']

PROVIDE_OPTIONS = []

def adjust_cfg(cfg):
    # Handle depending variables
    value = cfg.get_value('n')
    if type(value) == list:
        m = [1-1/n for n in value]
    else:
        m = 1 - 1/value
    cfg.set_value('m', value)

    # Set array indexes
    inner_points = cfg.get_value('inner_points')

    cfg.set_value('first_idx',    0)
    cfg.set_value('last_idx',     inner_points+1)
    cfg.set_value('mass_in_idx',  inner_points+2)
    cfg.set_value('s1_idx',       inner_points+3)
    cfg.set_value('s2_idx',       inner_points+4)
    cfg.set_value('mass_out_idx', inner_points+5)
    cfg.set_value('pq_idx',       inner_points+6)
    # total length of 'z' array (discretization points + s1,s2,mass_in,...)
    cfg.set_value('z_size',       inner_points+7)

    discretization_type = cfg.get_value('dtype')
    if discretization_type == 1:
        y = linspace(0, 1, inner_points + 2)
    else:
        raise NotImplementedError('For now only linear discretization is'
                                  ' implemented')

    cfg.set_value('y', y)
    cfg.set_value('y12', (y[1:]+y[:-1])/2.)

    dy = y[1:]-y[:-1]
    cfg.set_value('dy', dy)

    ldc1, ldc2, ldc3 = lagrangean_derivative_coefs(dy)
    cfg.set_value('ldc1', ldc1)
    cfg.set_value('ldc2', ldc2)
    cfg.set_value('ldc3', ldc3)

    cfg.set_value('calc_wm', cfg.get_value('calc_wm')
                             or cfg.get_value('calc_gc'))
