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

def load_modules_names(submodule):
    """
    Lists all experiment types defined in modules/<module_name>.info.py
    and returns a function for quick finding the module for given
    experiment type as input
    """

    # We keep the list of available modules in modules names, where
    # every modules_names[exp_type] points to given module name that
    # exports it. Modules itself are loaded only on demand.
    existing_modules = {}
    loaded_modules   = {}

    from os import listdir
    from sys import modules as sysmodules

    modules_names = listdir('modules')

    for module_name in modules_names:
        if module_name not in ['__pycache__', '__init__.py']:
            try:
                module_full_name = 'modules.' + module_name + '.info'
                __import__(module_full_name)

                info = sysmodules[module_full_name]

                for exp_id in info.types:
                    existing_modules[exp_id] = module_name
            except:
                print('Module load error:Module "%s" could not be loaded.'
                      ' Skipping.\n'  % module_name)

    def find_module(exp_type):
        if exp_type in existing_modules:
            module_name = existing_modules[exp_type]

            if not module_name in loaded_modules:
                module_full_name = 'modules.' + module_name + '.' + submodule
                __import__(module_full_name)
                loaded_modules[module_name] = sysmodules[module_full_name]

            return loaded_modules[module_name]
        else:
            print('Unknown experiment type: %s' % exp_type)
            print('Available experiments are:\n %s' % modules_names())
            raise ValueError()

    return find_module
