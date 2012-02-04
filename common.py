def read_model(cfg_filenames, base_cfg = {}, check_cfg_fn = None,
               adjust_cfg_fns = None):
    """
    Loads centrifuge inifiles (filenames - either a file or a list of files)
    and on the loaded data object calls the post_hook_fn (if set).
    """
    from config import read_configs, merge_cfgs, base_cfg
    from model_parameters import ModelParameters

    [cfg] = read_configs(inifilenames, preserve_sections_p=True, cfgs_merge=True)
    
    model_cfg = merge_cfgs(base_cfg, cfg)

    if adjust_cfg_fns:
        for fn in adjust_cfg_fns:
            fn(model_cfg)
    if check_cfg:
        check_cfg(model_cfg)

    model = ModelParameters(model_cfg)

    adjust_model_default(model, adjust_all = True)
    apply_functions(post_hook_fns, model)

    return model
