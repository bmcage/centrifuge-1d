def read_model(cfg_filenames, base_cfg = {}, check_cfg_fn = None,
               adjust_cfg_fns = None):
    """
    Loads centrifuge inifiles (filenames - either a file or a list of files)
    and on the loaded data object calls the post_hook_fn (if set).
    """
    from config import read_cfgs, merge_cfgs, base_cfg
    from model_parameters import ModelParameters

    [cfg] = read_cfgs(inifilenames, preserve_sections_p=True, cfgs_merge=True)
    
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

def load_modules_names(submodule)
    """
    Lists all experiment types defined in modules/<module_name>.info.py
    and returns a function for quick finding the module for given
    experiment type as input
    """

    # We keep the list of available modules in modules names, where
    # every modules_names[exp_type] points to given module name that
    # exports it. Modules itself are loaded only on demand.
    modules_names  = {}
    loaded_modules = {}

    import os.listdir as listdir

    modules_names = listdir
    for module_name in modules_names:
        if module_name not in ['__pycache__', '__init__']:
            try:
                info   = __import__('modules.' + module_name + '.info')

                for exp_id in info.types:
                    modules[exp_id] = module_name
            except:
                print('Module %s could not be loaded. Skipping.')

    def find_module(exp_type)
        if exp_type in modules_names:
            module_name = modules_names[exp_type]

            if not module_name in loaded_modules:
                module = __import__('modules.' + module_name + '.' + submodule)
                loaded_modules[module_name] = module

            return loaded_modules[module_name]
        else:
            print('Unknown experiment type: %s' % exp_type)
            print('Available experiments are:\n %s' % modules_names.keys())
            raise ValueError()
