from modules.shared.shared_functions import lagrangean_derivative_coefs
from numpy import linspace, power, asarray

PARENTAL_MODULES = ['base']

CONFIG_OPTIONS = {
        'mandatory' : ['inner_points', 'dtype', 'n', 'gamma', 'draw_graphs',
                       'h_init', 'porosity',
                       'calc_gc', 'calc_rm',
                       'rb_type'],
        'defaults'  : {},
        'dependent' : {'rb-2':
                         (lambda cfg: cfg.get_value('rb_type') == 2,
                          ['h_last']),
                        'rb-4':
                         (lambda cfg: cfg.get_value('rb_type') == 4,
                          ['h_last', 'dip_height'])},
        'optional'  : ['n1', 'gamma1', 'n2', 'gamma2'],
        'additional': ['m', 'y', 'y12', 'dy', 'ldc1', 'ldc2', 'ldc3',
                       'first_idx', 'last_idx', 's1_idx', 's2_idx',
                       'mass_in_idx', 'mass_out_idx', 'pq_idx', 'z_size',
                       'h_last']
        }

EXCLUDE_FROM_MODEL = ['dtype']

PROVIDE_OPTIONS = []

def adjust_cfg(cfg):
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

    if cfg.get_value('rb_type') == 4:
        omega2g = (power(asarray(cfg.get_value('omega'), dtype=float), 2)
                   /(2  * cfg.get_value('g')))
        rL = asarray(cfg.get_value('r0'))+asarray(cfg.get_value('l0'))
        h_last = omega2g * (power(rL
                                  + asarray(cfg.get_value('fl2'), dtype=float),
                                  2)
                            - power(rL - asarray(cfg.get_value('dip_height'),
                                                 dtype=float),
                                    2))
        cfg.set_value('h_last', h_last)
