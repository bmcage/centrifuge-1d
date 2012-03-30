PARENTAL_MODULES = ['base']

CONFIG_OPTIONS = {
        'mandatory' : ['inner_points', 'dtype', 'n', 'gamma', 'draw_graphs',
                       'h_init']
        'defaults'  : {}
        'dependent' : {}
        'optional'  : ['n1', 'gamma1', 'n2', 'gamma2']
        'additional': ['y', 'dy', 'ldc1', 'ldc2', 'ldc3',
                       'first_idx', 'last_idx', 's1_idx', 's2_idx',
                       'mass_in_idx', 'mass_out_idx', 'pq_idx', 'z_size',]
        }

CONFIG_ONLY_OPTIONS = []

IGNORE_OPTIONS = []
