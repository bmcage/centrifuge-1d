"""
This package implements access to configuration and defines
the ModelParameters class which stores the setting obtained
from the configuration files.
"""

try:
    import ConfigParser as configparser
except:
    import configparser
from numpy import pi
from os import listdir
from sys import modules as sysmodules
from math import sqrt

##################################################################
#                   ModulesManager class                         #
##################################################################

def get_ancestors(options_module):
    return options_module.PARENTAL_MODULES

class ModulesManager():
    def __init__(self):
        available_modules = listdir('modules')
        available_modules.remove('__pycache__')
        available_modules.remove('__init__.py')
        available_modules.remove('shared')

        loaded_modules = {}

        # module_names[exp_type]    -> module_name
        # module_names[module_name] -> module_name
        modules_names = {}

        for module_name in available_modules:
            try:
                submodule_info_name = 'modules.' + module_name + '.info'

                __import__(submodule_info_name)

                submodule_info = sysmodules[submodule_info_name]

                for exp_type in submodule_info.types:
                    modules_names[exp_type] = module_name

                modules_names[module_name] = module_name
            except:
                print('Module loading error:Submodule ''info'' of module "%s" '
                      'could not be loaded. Skipping.'  % module_name)

        self._loaded_modules = loaded_modules
        self._modules_names  = modules_names

    def traverse_ancestors(self, modname_or_exptype, apply_fn, submodule = '',
                           prehook = None, get_ancestors_fn = get_ancestors):
        """
          Recursively traverses the modules, preserving the order of ancestors
          it gets from 'get_ancestors(module)' function, and applies the
          function 'apply_fn(module)' on each (ancestor) module.
          If 'prehook(module)' is specified, this function is called before
          going into recursion.
        """

        def traverse(module):
            ancestors_list = get_ancestors_fn(module)

            for ancestor_name in ancestors_list:
                ancestor_module = self.find_module(ancestor_name, submodule)

                if prehook:
                    prehook(ancestor_module)

                if not traverse(ancestor_module):
                    return False

            return apply_fn(module)

        module = self.find_module(modname_or_exptype, submodule)

        return traverse(module)

    def find_module(self, modname_or_exptype, submodule=''):
        """
          Return an module determined by it's name or the type of experiment.
          Loading of modules is done only on demand.
        """

        if not modname_or_exptype in self._modules_names:
            print('\n\n',modname_or_exptype, self._modules_names,'\n\n')
            print('\nFind module error: Unknown module name or experiment '
                  'type: "%s".\nAvailable modules are: %s'
                  % (modname_or_exptype, list(self._modules_names.keys())))
            raise
            exit(1)

        if submodule:
            module_name = (self._modules_names[modname_or_exptype]
                           + '.' + submodule)
        else:
             module_name = self._modules_names[modname_or_exptype]

        if not module_name in self._loaded_modules:
            module_full_name = 'modules.' + module_name
            __import__(module_full_name)

            module = sysmodules[module_full_name]
            self._loaded_modules[module_name] = module
        else:
            module = self._loaded_modules[module_name]

        return module

##################################################################
#                   Configuration class                          #
##################################################################

def flatten(cfg):
    cfg_dict = cfg._cfg_dict

    def flatten_dict(base_dict, dict_tree):
        for (key, value) in dict_tree.items():
            if type(value) == dict:
                dict_flatten(base_dict, value)
            else:
                base_dict[key] = value

    flattened_cfg = Configuration()
    flattened_cfg._cfg_dict = flatten_dict({}, cfg)

    return flattened_cfg

def parse_value(str_value):
    """
      Given a value as string, tries to convert to it's correspending type.
      May be called recursively in the case of nested structures.
    """
    try:
        raw_value = str_value.strip()

        if raw_value[0] == "[" and raw_value[-1] == "]":
            return [parse_value(item) for item in raw_value[1:-1].split(",")]
        elif ((raw_value[0] == "'" and raw_value[-1] == "'")
              or (raw_value[0] == '"' and raw_value[-1] == '"')):
            return raw_value[1:-1]
        elif raw_value == 'True':
            return True
        elif raw_value == 'False':
            return False
        elif "." in raw_value or "e" in raw_value or "E" in raw_value:
            return float(raw_value)
        elif raw_value[0] == "(" and raw_value[-1] == ")":
            return \
              tuple([parse_value(item) for item in raw_value[1:-1].split(",")])
        else:
            return int(raw_value)

    except:
        print('Error:Could not parse value: ', raw_value, '\nExiting...')
        exit(1)

class Configuration:
    def __init__(self, preserve_sections_p = False):

        self._preserve_sections_p = preserve_sections_p
        self._cfg_dict = {}

    def set_value(self, key, value, section = None):
        cfg_dict = self._cfg_dict

        if self._preserve_sections_p:
            if not section:
                print('set_value error: Section not specified: ', section)
                exit(1)
            cfg_dict[section][key] = value
        else:
            cfg_dict[key] = value

    def get_value(self, key, section = None):

        cfg_dict = self._cfg_dict

        if self._preserve_sections_p and (not section
                                          or not section in cfg_dict):
            print('get_value error: Section not present in config: ', section)
            exit(1)

        if self._preserve_sections_p:
            if key in cfg_dict[section]:
                return cfg_dict[section][key]
            else:
                return None
        else:
            if key in cfg_dict:
                return cfg_dict[key]
            else:
                return None

    def missing_options(self, options, section = None):

        if self._preserve_sections_p:
            if (not section) or (not section in self._cfg_dict):
                print('missing_options error: Section not present in config: ',
                      section)
                exit(1)

            cfg_dict = self._cfg_dict[section]
        else:
            cfg_dict = self._cfg_dict

        result = []

        for option in options:
            if not option in cfg_dict:
                result.append(option)

        return result

    def list_options(self):
        """
          Return a list of all options specified in config.
        """
        return list(self._cfg_dict.keys())

    def echo(self):
        def echo_section(section_dict):
            for (name, value) in sorted(section_dict.items()):
                print('%-12s = %s' % (name, value))

        cfg_dict = self._cfg_dict

        if self._preserve_sections_p:
            for (section, section_value) in cfg_dict.items():
                print('\n[%s]', section)
                echo_section(section_value)
        else:
            print()
            echo_section(cfg_dict)

    def set_defaults(self, modman):
        """
          Set options that are missing in current configuration, but have
          specified a default value in modules.
        """

        cfg_dict = self._cfg_dict

        def set_defaults_from_module(options_module):
            defaults_options = options_module.CONFIG_OPTIONS['defaults']

            if not self._preserve_sections_p:
                options = {'foo': defaults_options}

            for (section, section_content) in options.items():
                for (option, value) in section_content.items():
                    if not option in cfg_dict:
                        self.set_value(option, value, section)

        exp_type = self.get_value('exp_type')

        modman.traverse_ancestors(exp_type, set_defaults_from_module,
                                  submodule='options')

        return self

    def adjust_cfg(self, modman):

        def module_adjust_cfg(module):
            if hasattr(module, 'adjust_cfg'):
                module.adjust_cfg(self)

        exp_type = self.get_value('exp_type')

        modman.traverse_ancestors(exp_type, module_adjust_cfg,
                                  submodule='options')

    def merge(self, *cfgs):
        """
          Merge all cfgs into current configuration; if the same values appears
          in more than one of the cfgs, the one from the last cfg applies.
        """
        cfg_dict = self._cfg_dict
        preserve_sections_p = self._preserve_sections_p

        for cfg in cfgs:
            if not cfg: continue

            if not (preserve_sections_p == cfg._preserve_sections_p):
                print('MergeCFG Error: Configurations do have same'
                      ' ''preserve_sections_p'' flag. Cannot merge.')
                exit(1)

            if preserve_sections_p:
                for (section, section_value) in cfg.items():
                    if not section in cfg_dict:
                        cfg_dict[section] = {}

                    cfg_dict[section].update(section_value)
            else:
                cfg_dict.update(cfg._cfg_dict)

        return self

    def read_from_files(self, *cfgs_filenames):
        """
          Reads configuration from .ini files 'cfgs_filenames'.
          If preserve_sections_p is false, removes sections and leaves only
          values from all sections; if two sections contained
          the same values, the latter will be kept.
        """

        # read config files
        parser     = configparser.ConfigParser()
        read_files = parser.read(cfgs_filenames)

        if (len(read_files) != len(cfgs_filenames)):
            print('Warning: from expected files: ', str(cfgs_filenames),
                  'were successfully parsed only: ', str(read_files))

        # Write data from parser to configuration
        cfg_dict = self._cfg_dict
        preserve_sections_p = self._preserve_sections_p

        for psection in parser.sections():
            section = psection.lower()

            if preserve_sections_p and not (section in cfg_dict):
                cfg_dict[section] = {}

            for option in parser.options(psection):
                raw_value = parser.get(psection, option).strip()

                value = parse_value(raw_value)

                if preserve_sections_p:
                    cfg_dict[section][option]=value
                else:
                    cfg_dict[option]=value

        return self

    def is_valid(self, modman):

        ignored_options = []
        # TODO: add ignore_options
        #ignored_options = set(ignore_options)


        def check_options(options):
            required_options = set(options)
            required_options.difference_update(ignored_options)

            missing_options = self.missing_options(required_options)
            if missing_options:
                print('Following required options are not present: ')
                for option in missing_options:
                    print('  ', option)

                    return False
            return True

        alien_options = set(self.list_options())

        def check_cfg(options_module):
            config_parameters = options_module.CONFIG_OPTIONS

            if not check_options(config_parameters['mandatory']):
                return False

            # remove known options from the list
            alien_options.difference_update(set(config_parameters['mandatory']))

            if 'dependent' in config_parameters:
                dependent_options = config_parameters['dependent']

                for (test_fn, dependent_options) \
                  in config_parameters['dependent'].values():

                    if test_fn(self):
                        if not check_options(self, dependent_options):
                            return False
                    alien_options.difference_update(set(dependent_options))

            if 'optional' in config_parameters:
                optional_options = config_parameters['optional']
                alien_options.difference_update(set(optional_options))

            if 'defaults' in config_parameters:
                defaults_options = config_parameters['defaults']
                alien_options.difference_update(set(defaults_options))

            return True

        exp_type = self.get_value('exp_type')

        if not exp_type:
            print('Missing option: configuration does not specify ''exp_type'' '
                  'option. Cannot validate configuration.\nCurrent '
                  'configuration is:')
            self.echo()
            print('\nExiting...')
            exit(1)

        if not modman.traverse_ancestors(exp_type, check_cfg,
                                         submodule='options'):
            print('ddddd')
            return False

        if alien_options:
            print('\nFound following alien options in configuration:')
            for option in alien_options:
                print('  ', option)
            return False

        return True

##################################################################
#                 ModelParameters class                          #
##################################################################

def rpm2radps(x):
    """
      Converts rpm to rad.s^{-1}
    """
    # rpm->rad.s-1:  omega_radps = (2pi)*omega_rps/60
    return x * pi/ 30.0

class ModelParameters:
    """
    Parameters class for the centrifuge simulation.
    """

    def __init__(self, cfg, modman):
        self._cfg = cfg
        self._iterable_parameters = {}
        self.iteration            = 0
        self.iterations           = 0

        self.iterations = len(cfg.get_value('duration'))

        exclude_options = []

        def update_exclude_option(options_module):
            #global exclude_options
            exclude_options.extend(options_module.EXCLUDE_FROM_MODEL)

        def set_options(options_list):
             for option in options_list:
                if not option in exclude_options:
                    self.set_value(option, cfg.get_value(option))

        def set_options_from_config(options_module):
            config_options = options_module.CONFIG_OPTIONS

            set_options(config_options['mandatory'])
            set_options(config_options['defaults'].keys())
            set_options(config_options['additional'])

            optional_options = set(config_options['optional'])
            cfg_missing_options  = set(cfg.missing_options(optional_options))
            optional_options.difference_update(cfg_missing_options)
            set_options(list(optional_options))

            for (test_fn, dependent_options) \
              in config_options['dependent'].values():

                if test_fn(cfg):
                    set_options(dependent_options)

        modman.traverse_ancestors(cfg.get_value('exp_type'),
                                  set_options_from_config,
                                  prehook=update_exclude_option,
                                  submodule='options')

    def set_value(self, key, value):
        """
          Set the value of parameter given as 'key'.
          Performs the check of the 'value' and conversion to correct units
          when needed. Updates also the dependent variables.
        """

        if key in ['omega_fall', 'm']:
            raise ValueError('Parameter %s is not allowed to set directly.'
                             % key)

        # Keep self._itarable_parameters up-to-date; if we set a list-type value
        # should be stored, if an atom, should be removed
        if (type(value) == list):
            # Handle specialy-treated variables, for which the value
            # needs to be somehow trasformed
            if key in ['omega', 'omega_start', 'omega_end']:
                self._iterable_parameters[key] = \
                  [rpm2radps(omega) for omega in value]
                return

            # Ok, no special treatment, just set it then
            self._iterable_parameters[key] = value

            # Update depending variables (value of which depends on the value
            # of variables in configuration)
            if key == 'n':
                self._iterable_parameters['m'] = [1-1/n for n in value]
            if key in ['r0_fall', 'g']:

                if key == 'r0_fall':
                    r0_fall_list = value

                    if 'g' in self._iterable_parameters:
                        g_list = self._iterable_parameters['g']
                    elif hasattr(self, 'g'):
                        g_list = [self.g for r0_fall in r0_fall_list]
                    else:
                        return
                else:
                    g_list = value

                    if 'r0_fall' in self._iterable_parameters:
                        r0_fall_list = self._iterable_parameters['r0_fall']
                    elif hasattr(self, 'r0_fall'):
                        r0_fall_list = [self.r0_fall for g in g_list]
                    else:
                        return

                omega_fall = [sqrt(g/r0_fall)
                              for (g, r0_fall) in zip(g_list, r0_fall_list)]
                self._iterable_parameters['omega_fall'] = omega_fall
        else:
            # Now is value an atom, so remove it from iterables if present
            if key in self._iterable_parameters:
                del(self._iterable_parameters[key])

            # Handle variables that need to transform the supplied value
            if key in ['omega', 'omega_start', 'omega_end']:
                setattr(self, key, rpm2radps(value))
                return

            # Handle the rest of supplied variable
            setattr(self, key, value)

            # Initialize depending variables
            if key == 'n':
                setattr(self, 'm', m = 1-1/value)
            elif key in ['r0_fall', 'g']:
                if hasattr(self, 'g') and hasattr(self, 'r0_fall'):
                    setattr(self, 'omega_fall', sqrt(self.g/self.r0_fall))

    def next_iteration(self):
        """
          Assign the next value of the parameters that were given as type list
        """
        i = self.iteration

        cfg = self._cfg

        for (key, value) in self._iterable_parameters.items():
            setattr(self, key, value[i])

        self.iteration = i+1

        return (i < self.iterations)

    def echo(self, iterable_only=False):
        """
          Print the parameters stored in the model.
          If 'iterable_only' is True, prints only the current values of the
          variables that are iterated by the 'next_iteration()' function.
        """

        iterable_options = self._iterable_parameters
        print('\nIterable parameters:')
        for option in sorted(iterable_options):
            print('  %-12s = %s' % (option, iterable_options[option]))

        if not iterable_only:
            options = dict(vars(self))
            for internal_opt in ['_cfg', '_iterable_parameters', 'iteration',
                                 'iterations']:
                del(options[internal_opt])

            print('\nConstant parameters:')
            for option in sorted(options):
                print('  %-12s = %s' % (option, getattr(self, option)))

            print('\nParameters set internally:')
            options = ['iteration', 'iterations']
            for option in options:
                print('  %-12s = %s' % (option, getattr(self, option)))
