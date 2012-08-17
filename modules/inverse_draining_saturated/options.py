PARENTAL_MODULES = ['direct_draining_saturated']

CONFIG_OPTIONS = {
        'mandatory' : ['inv_init_params'],
        'defaults'  : {'optimfn': 'leastsq'},
        'dependent' : {'inv_params':
                       (lambda cfg: len(cfg.get_value('inv_init_params')) == 2,
                        ['ks'])},
        'optional'  : ['wl_out', 'gc1', 'rm1', 'inv_ubounds', 'inv_lbounds'],
        'additional': ['calc_wl_out', 'params_ubounds', 'params_lbounds',
                       'ks_inv_scale']
        }

EXCLUDE_FROM_MODEL = ['inv_ubounds', 'inv_lbounds']

PROVIDE_OPTIONS = [(lambda cfg: len(cfg.get_value('inv_init_params')) == 3, ['ks']),
                   'n', 'gamma', 'calc_gc', 'calc_rm']

NONITERABLE_LIST_OPTIONS = ['inv_init_params']

def check_cfg(cfg):
    inv_params_len = len(cfg.get_value('inv_init_params'))
    if not ((inv_params_len == 3) or (inv_params_len == 2)):
        print('Option ''inv_params'' should be of length 3 or 2.')
        return False

    inv_ubounds = cfg.get_value('inv_ubounds')
    inv_lbounds = cfg.get_value('inv_lbounds')
    if (((not inv_ubounds  is None) and (len(inv_ubounds) != inv_params_len))
        or ((not inv_lbounds  is None)
            and (len(inv_lbounds) != inv_params_len))):
        print('The length of options ''inv_lbounds'' and ''inv_lbounds'' must '
              'be the same as the length of ''inv_params'' option.')
        return False

    if ((cfg.get_value('wl_out') is None) and (cfg.get_value('gc1') is None)
        and (cfg.get_value('rm1') is None)):
        print('No measured data (wl_out, gc1, rm1) is specified. Cannot run '
              'inverse problem')
        return False

    return True

def adjust_cfg(cfg):
    from numpy import inf, asarray, power, trunc, log10

    cfg.set_value('calc_gc', not cfg.get_value('gc1') is None)
    cfg.set_value('calc_rm', not cfg.get_value('rm1') is None)
    cfg.set_value('calc_wl_out', not cfg.get_value('wl_out') is None)

    params_ubounds = cfg.get_value('inv_ubounds')
    params_lbounds = cfg.get_value('inv_lbounds')

    determine_all = (len(cfg.get_value('inv_init_params')) == 3)
    if determine_all:
        params_names = ['ks', 'n', 'gamma']
        if params_ubounds is None:
            params_ubounds = asarray([inf, inf, 0.])
        if params_lbounds is None:
            params_lbounds = asarray([0., 1., -inf])

        ks_init = cfg.get_value('inv_init_params')[0]
        cfg.set_value('ks_inv_scale', power(10, -trunc(log10(ks_init))))

    else:
        params_names = ['n', 'gamma']
        if params_ubounds is None:
            params_ubounds = asarray([inf, 0.])
        if params_lbounds is None:
            params_lbounds = asarray([1., -inf])
        cfg.set_value('ks_inv_scale', 1.0)

    cfg.set_value('params_ubounds', dict(zip(params_names, params_ubounds)))
    cfg.set_value('params_lbounds', dict(zip(params_names, params_lbounds)))
