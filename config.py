"""
This package implements access to configuration and defines
the ModelParameters class which stores the setting obtained
from the configuration files.
"""

try:
    import ConfigParser as configparser
except:
    import configparser
from numpy import pi, inf
from os import listdir
from sys import modules as sysmodules
from math import sqrt

##################################################################
#                   ModulesManager class                         #
##################################################################

def get_ancestors(options_module):
    return options_module.PARENTAL_MODULES


def rpm2radps(x):
    """
      Converts rpm to rad.s^{-1}
    """
    # rpm->rad.s-1:  omega_radps = (2pi)*omega_rps/60
    return x * pi/ 30.0

def calc_omega_fall(r0_fall, g):
    is_list_r0_fall = type(r0_fall) == list
    is_list_g       = type(g) == list

    if is_list_r0_fall and is_list_g:
        omega_fall = [sqrt(g_f/r0_f) for (g_f, r0_f) in zip(g, r0_fall)]
    elif is_list_r0_fall:
        omega_fall = [sqrt(g/r0_f) for r0_f in r0_fall]
    elif is_list_g:
        omega_fall = [sqrt(g_f/r0_fall) for g_f in g]
    else:
        omega_fall = sqrt(g/r0_fall)

    return omega_fall

class ModulesManager():
    def __init__(self):
        available_modules = listdir('modules')
        available_modules.remove('__pycache__')
        available_modules.remove('__init__.py')

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

        self._available_modules = available_modules
        self._loaded_modules = loaded_modules
        self._modules_names  = modules_names

    def echo(self, verbose=False):
        """
          Display information about available modules and experiment types
          provided by every module.
          If 'verbose' is set to 'True', display also module description.
        """
        find_module = self.find_module

        print('\n\nCentrifuge modules:')

        for module_name in sorted(self._available_modules):
            module = find_module(module_name, submodule='info')
            print('\n  %s:' % module_name)
            print('        Experiment types:', module.types)
            if verbose:
                print('        Description:')
                print('            ', module.desc)

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

            if prehook:
                prehook(module)

            ancestors_list = get_ancestors_fn(module)

            for ancestor_name in ancestors_list:
                ancestor_module = self.find_module(ancestor_name, submodule)

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
            value = []
            for item in raw_value[1:-1].split(","):
                one_value = parse_value(item)
                if type(one_value) == list:
                    value.extend(one_value)
                else:
                    value.append(one_value)
            return value
        elif ((raw_value[0] == "'" and raw_value[-1] == "'")
              or (raw_value[0] == '"' and raw_value[-1] == '"')):
            return raw_value[1:-1]
        elif raw_value == 'True':
            return True
        elif raw_value == 'False':
            return False
        elif raw_value == 'None':
            return None
        elif raw_value == 'inf':
            return inf
        elif raw_value == '-inf':
            return -inf
        elif '*' in raw_value:
            [raw_mul, raw_val] = raw_value.split("*")
            mul = parse_value(raw_mul)
            val = parse_value(raw_val)
            return [val for i in range(mul)]
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

        if self._preserve_sections_p:
            if not section:
                print('set_value error: Section not specified: ', section)
                exit(1)
            cfg_dict = self._cfg_dict[section]
        else:
            cfg_dict = self._cfg_dict

        cfg_dict[key] = value

        # Handle depending variables
        if key == 'n':
            if type(value) == list:
                m = [1-1/n for n in value]
            else:
                m = 1 - 1/value
            cfg_dict['m'] = m

        elif key in ['r0_fall', 'g']:
            if ((key == 'r0_fall' and 'g' in cfg_dict)
                or (key == 'g' and 'r0_fall' in cfg_dict)):

                cfg_dict['omega_fall'] = \
                  calc_omega_fall(cfg_dict['r0_fall'], cfg_dict['g'])

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

        blacklisted_options = set([])

        def update_special_options(options_module):
            if hasattr(options_module, 'BLACKLIST_OPTIONS'):
                blacklisted_options.update(options_module.BLACKLIST_OPTIONS)

        def set_defaults_from_module(options_module):
            defaults_options = options_module.CONFIG_OPTIONS['defaults']

            if not self._preserve_sections_p:
                options = {'foo': defaults_options}

            for (section, section_content) in options.items():
                for (option, value) in section_content.items():
                    if ((not option in cfg_dict)
                        and (not option in blacklisted_options)):
                        self.set_value(option, value, section)

            return True

        exp_type = self.get_value('exp_type')

        modman.traverse_ancestors(exp_type, set_defaults_from_module,
                                  submodule='options',
                                  prehook=update_special_options)

        return self

    def adjust_cfg(self, modman):

        def module_adjust_cfg(module):
            if hasattr(module, 'adjust_cfg'):
                module.adjust_cfg(self)
                return True

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

        def transform_value(key, value):
            """
              Handle specialy-treated values, that need to be transformed
              (e.g. into correct units)
            """
            if key in ['omega', 'omega_start', 'omega_end']:
                if type(value) == list:
                  return [rpm2radps(omega) for omega in value]
                else:
                    return rpm2radps(value)
            else:
                return value

        # read config files
        parser     = configparser.ConfigParser()
        try:
            read_files = parser.read(cfgs_filenames)
        except configparser.DuplicateOptionError as E:
            print(E)
            exit(0)

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

                value = transform_value(option, parse_value(raw_value))
                self.set_value(option, value, section)

        return self

    def is_valid(self, modman):

        provided_options = set([])
        blacklisted_options = set([])

        def update_special_options(options_module):
            if hasattr(options_module, 'PROVIDE_OPTIONS'):
                provided_options.update(options_module.PROVIDE_OPTIONS)
            if hasattr(options_module, 'BLACKLIST_OPTIONS'):
                blacklisted_options.update(options_module.BLACKLIST_OPTIONS)

        alien_options = set(self.list_options())

        def handle_special_options(options):
            """
              From 'options' remove provided and blacklisted options and
              modify accordingly the 'alien_options' list.
            """

            normal_options = set(options)
            normal_options.difference_update(provided_options)
            normal_options.difference_update(blacklisted_options)

            alien_options.difference_update(normal_options)

            return normal_options

        def check_options(options):
            """
              Check whether all supplied 'options' in present configuration.
              current configuration.
            """
            missing_options = self.missing_options(options)
            if missing_options:
                print('Following required options are not present: ')
                for option in missing_options:
                    print('  ', option)

                return False
            return True

        def primary_cfg_check(options_module):
            config_parameters = options_module.CONFIG_OPTIONS

            required_options = \
              handle_special_options(config_parameters['mandatory'])

            if not check_options(required_options):
                return False

            if 'dependent' in config_parameters:
                for (test_fn, dependent_options) \
                  in config_parameters['dependent'].values():

                    if test_fn(self):
                        required_options = \
                          handle_special_options(dependent_options)
                        if not check_options(required_options):
                            return False

            if 'optional' in config_parameters:
                handle_special_options(config_parameters['optional'])

            if 'defaults' in config_parameters:
                handle_special_options(config_parameters['defaults'])

            if 'additional' in config_parameters:
                handle_special_options(config_parameters['additional'])

            return True

        def custom_cfg_check(options_module):
            if hasattr(options_module, 'check_cfg'):
                return options_module.check_cfg(self)
            else:
                return True


        exp_type = self.get_value('exp_type')

        if not exp_type:
            print('Missing option: configuration does not specify ''exp_type'' '
                  'option. Cannot validate configuration.\nCurrent '
                  'configuration is:')
            self.echo()
            print('\nExiting...')
            exit(1)

        if not modman.traverse_ancestors(exp_type, primary_cfg_check,
                                         submodule='options',
                                         prehook=update_special_options):
            return False

        if not modman.traverse_ancestors(exp_type, custom_cfg_check,
                                         submodule='options'):
            return False

        if alien_options:
            print('\nFound following alien options in configuration:')

            provided_aliens = alien_options.intersection(provided_options)
            if not provided_aliens.issubset([]):
                print('\n  Options found in configuration, but PROVIDED by '
                      'a centrifuge module:')
                for option in provided_aliens:
                    print('    ', option)

                alien_options.difference_update(provided_aliens)

            blacklisted_aliens = alien_options.intersection(blacklisted_options)
            if not blacklisted_aliens.issubset([]):
                print('\n  Options found in configuration, but BLACKLISTED by '
                      'a centrifuge module:')
                for option in blacklisted_aliens:
                    print('    ', option)

                    alien_options.difference_update(blacklisted_aliens)

            if alien_options:
                print('\n Options found in configuration, but not specified'
                      'by a module:')
                for option in alien_options:
                    print('  ', option)

            return False

        return True

##################################################################
#                 ModelParameters class                          #
##################################################################

class ModelParameters:
    """
    Parameters class for the centrifuge simulation.
    """

    def __init__(self, cfg, modman):
        self._cfg = cfg
        self._iterable_parameters = {}
        self.iteration            = 0
        self.iterations           = 0
        self._atomic_options      = []

        exclude_options = []
        atoms           = self._atomic_options

        def update_options_list(options_module):
            #global exclude_options
            if hasattr(options_module, 'EXCLUDE_FROM_MODEL'):
                exclude_options.extend(options_module.EXCLUDE_FROM_MODEL)

            if hasattr(options_module, 'NONITERABLE_LIST_OPTIONS'):
                atoms.extend(options_module.NONITERABLE_LIST_OPTIONS)

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

            return True

        modman.traverse_ancestors(cfg.get_value('exp_type'),
                                  set_options_from_config,
                                  prehook=update_options_list,
                                  submodule='options')

        if self._iterable_parameters:
            iterations = -1

            for (key, value) in self._iterable_parameters.items():
                if iterations == -1:
                    ref_key = key
                    iterations = len(value)
                else:
                    if len(value) != iterations:
                        print('Option list ''%s'' has length: %i\n'
                              'Reference list ''%s'' has length: %i\n'
                              'Cannot iterate over lists of unequal length.'
                              'Aborting...'
                              % (key, len(value), ref_key, iterations))
                        exit(1)
        else:
            iterations = 1

        self.iterations = iterations

    def get_iterable_value(self, key):
        if key in self._iterable_parameters:
            return self._iterable_parameters[key]
        else:
             return None

    def set_value(self, key, value):
        """
          Set the value of parameter given as 'key'.
        """

        # Keep self._itarable_parameters up-to-date; if we set a list-type value
        # should be stored, if an atom, should be removed
        if (type(value) == list) and (not key in self._atomic_options):
            self._iterable_parameters[key] = value
        else:
            # Now is value an atom, so remove it from iterables if present
            if key in self._iterable_parameters:
                del(self._iterable_parameters[key])

            # Handle the rest of supplied variable
            setattr(self, key, value)

    def next_iteration(self):
        """
          Assign the next value of the parameters that were given as type list
        """
        i = self.iteration
        if i == self.iterations: return False

        cfg = self._cfg

        for (key, value) in self._iterable_parameters.items():
            setattr(self, key, value[i])

        self.iteration = i+1

        return True

    def init_iteration(self):
        self.iteration = 0

        for (key, value) in self._iterable_parameters.items():
            setattr(self, key, value[0])

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
