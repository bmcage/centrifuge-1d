"""
This package implements access to configuration and defines
the ModelParameters class which stores the setting obtained
from the configuration files.
"""

try:
    import ConfigParser as configparser
except:
    import configparser
from numpy import inf, concatenate, cumsum, asarray, any as np_any
from modules.shared.functions import has_data
from os import listdir
from sys import modules as sysmodules
from types import MethodType

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
        elif raw_value[0] == "(" and raw_value[-1] == ")":
            return \
              tuple([parse_value(item) for item in raw_value[1:-1].split(",")])
        elif raw_value[0] == "{" and raw_value[-1] == "}":
            # ...what an interesting parsing... but who is going to implement
            # this again?
            value = eval(raw_value)
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
        else:
            return int(raw_value)

    except:
        print('Error:Could not parse value: ', raw_value, '\nExiting...')
        exit(1)

class Configuration:
    def __init__(self, preserve_sections_p = False):

        self._preserve_sections_p = preserve_sections_p
        self._cfg_dict = {}
        self._config_definition = None

    def set_value(self, key, value, section = None):

        if self._preserve_sections_p:
            if not section:
                print('set_value error: Section not specified: ', section)
                exit(1)
            cfg_dict = self._cfg_dict[section]
        else:
            cfg_dict = self._cfg_dict

        cfg_dict[key] = value

    def get_value(self, key, section = None, not_found=None):

        cfg_dict = self._cfg_dict

        if self._preserve_sections_p and (not section
                                          or not section in cfg_dict):
            print('get_value error: Section not present in config: ', section)
            exit(1)

        if self._preserve_sections_p:
            if key in cfg_dict[section]:
                return cfg_dict[section][key]
            else:
                return not_found
        else:
            if key in cfg_dict:
                return cfg_dict[key]
            else:
                return not_found

    def del_value(self, key, section = None):
        if self._preserve_sections_p:
            if not section:
                print('set_value error: Section not specified: ', section)
                exit(1)
            cfg_dict = self._cfg_dict[section]
        else:
            cfg_dict = self._cfg_dict

        del cfg_dict[key]

    def iterate_values(self, section = None):
        if self._preserve_sections_p and section:
            print('cfg:iterate_values: preserving sections is not implemented.')
            exit(1)

        return map(lambda value: value, self._cfg_dict.items())

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

        if not self._config_definition:
            self._load_config_definition(modman)

        defaults_options    = self._config_definition['default']
        blacklisted_options = self._config_definition['blacklisted']

        if not self._preserve_sections_p:
            section = 'foo'
        else:
            raise NotImplementedError('Section preserving is not implemented.')

        for (option, value) in defaults_options:
            if ((not option in cfg_dict)
                and (not option in blacklisted_options)):
                self.set_value(option, value, section)

        return self

    def _group(self):
        def _measurements_times(cfg):
            t_duration    = cfg.get_value('duration', not_found=None)
            t_fh_duration = cfg.get_value('fh_duration', not_found=None)
            dec_duration  = cfg.get_value('deceleration_duration',
                                          not_found=None)

            if t_duration is None:
                t_duration = 0.0
            else:
                t_duration = asarray(t_duration, dtype=float)
            if dec_duration  is None: dec_duration = 0.0
            if cfg.get_value('include_acceleration'):
                t_duration[:] = t_duration + dec_duration

            if t_fh_duration is None:
                t_fh_duration = 0.0
            else:
                t_fh_duration = asarray(t_fh_duration, dtype=float)

            duration_times = t_duration + t_fh_duration
            if not np_any(duration_times): return None # no times were specified

            stop_times = cumsum(concatenate(([0.0], duration_times)))

            return stop_times

        def _group_measurments(cfg):
            measurements_names = ('MI', 'MO', 'GC', 'RM', 'theta')

            measurements = {}
            measurements_weights = {}

            t = _measurements_times(cfg)
            if not t is None:
                measurements['t'] = t
                t_meas = t[1:]

            for name in measurements_names:
                if name == 'MI':
                    iname = 'wl1'
                elif name == 'MO':
                    iname = 'wl_out'
                else:
                    iname = name

                value = cfg.get_value(iname, not_found=None)
                if value == None: continue

                cfg.del_value(iname)
                if name == 'MO':
                    value = cumsum(asarray(value, dtype=float))

                if name == 'theta':
                    measurements[name] = (value, cfg.get_value('p'))
                else:
                    measurements[name] = (t_meas, value)

                iname_w = iname + '_weights'
                value_w = cfg.get_value(iname_w, not_found=None)
                if value_w == None: continue

                self.del_value(iname_w)
                measurements_weights[name] = value_w

            cfg.set_value('measurements', measurements)
            cfg.set_value('measurements_weights', measurements_weights)

        # Function body
        _group_measurments(self)

    def adjust_cfg(self, modman):

        def module_preadjust_cfg(module):
            if hasattr(module, 'prior_adjust_cfg'):
                module.prior_adjust_cfg(self)
            return True

        def module_adjust_cfg(module):
            if hasattr(module, 'adjust_cfg'):
                module.adjust_cfg(self)
            return True

        exp_type = self.get_value('exp_type')

        modman.traverse_ancestors(exp_type, module_adjust_cfg,
                                  prehook=module_preadjust_cfg,
                                  submodule='options')

        self._group()

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
            return value

        cfg_dict = self._cfg_dict
        preserve_sections_p = self._preserve_sections_p

        for fname in cfgs_filenames:
            # read config files
            parser   = configparser.ConfigParser()
            try:
                read_files = parser.read(fname)
            except configparser.DuplicateOptionError as E:
                print(E)
                exit(0)

            # Write data from parser to configuration
            for psection in parser.sections():
                section = psection.lower()

                if preserve_sections_p and not (section in cfg_dict):
                    cfg_dict[section] = {}

                for option in parser.options(psection):
                    raw_value = parser.get(psection, option).strip()

                    value = transform_value(option, parse_value(raw_value))
                    self.set_value(option, value, section)

        # we loaded new configuration so discard previously loaded definition
        self._config_definition = None

        return self

    def _load_config_definition(self, modman):

        def parse_options(options):
            required_options = []
            default_options = []

            def classify_options(options):
                for option in options:
                    option_type = type(option)
                    if option_type == str:
                        required_options.append(option)
                    elif (option_type == list) or (option_type == tuple):
                        if hasattr(option[0], '__call__'):
                            if option[0](self):
                                classify_options(option[1])
                            elif len(option) == 3:
                                classify_options(option[2])
                        else:
                            default_options.append(option)
                    elif hasattr(option, '__call__'):
                        classify_options(option(self))
                    else:
                        print('Unknown option type:', option)
                        exit(1)

            classify_options(options)

            return required_options, default_options

        required_options    = set([])
        default_options     = []
        provided_options    = set([])
        blacklisted_options = set([])
        internal_options    = set([])
        excluded_options    = set([])
        iterable_options    = set([])

        def read_options_definition(options_module):
            if hasattr(options_module, 'CONFIG_OPTIONS'):
                (required, defaults) = parse_options(options_module.CONFIG_OPTIONS)
                required_options.update(required)
                default_options.extend(defaults)

            options_name = 'PROVIDE_OPTIONS'
            if hasattr(options_module, options_name):
                (options, defaults) = \
                  parse_options(getattr(options_module, options_name))
                if defaults:
                    raise ValueError("Options in '%s' cannot have a "
                                     "default value: " % options_name,
                                     getattr(options_module, options_name))
                provided_options.update(options)

            options_name = 'BLACKLIST_OPTIONS'
            if hasattr(options_module, options_name):
                (options, defaults) = \
                  parse_options(getattr(options_module, options_name))
                if defaults:
                    raise ValueError("Options in '%s' cannot have a "
                                     "default value: " % options_name,
                                     getattr(options_module, options_name))
                blacklisted_options.update(options)

            # options_name = 'INTERNAL_OPTIONS'
            # if hasattr(options_module, options_name):
            #     (options, defaults) = \
            #       parse_options(getattr(options_module, options_name))
            #     if defaults:
            #         raise ValueError("Options in '%s' cannot have a "
            #                          "default value: " % options_name,
            #                          getattr(options_module, options_name))
            #     internal_options.update(options)


            options_name = 'EXCLUDE_FROM_MODEL'
            if hasattr(options_module, options_name):
                (options, defaults) = \
                  parse_options(getattr(options_module, options_name))
                if defaults:
                    raise ValueError("Options in '%s' cannot have a "
                                     "default value: " % options_name,
                                     getattr(options_module, options_name))
                excluded_options.update(options)

            options_name = 'OPTIONS_ITERABLE_LISTS'
            if hasattr(options_module, options_name):
                (options, defaults) = \
                  parse_options(getattr(options_module, options_name))
                if defaults:
                    raise ValueError("Options in '%s' cannot have a "
                                     'default value: ' % options_name,
                                     getattr(options_module, options_name))
                iterable_options.update(options)

            return True

        exp_type = self.get_value('exp_type')

        if not exp_type:
            print('Missing option: configuration does not specify ''exp_type'' '
                  'option. Cannot validate configuration.\nCurrent '
                  'configuration is:')
            self.echo()
            print('\nExiting...')
            exit(1)

        if not modman.traverse_ancestors(exp_type, read_options_definition,
                                         submodule='options'):
            print('Could not parse options definitions... Exiting...')
            exit(1)

        # from required options remove blacklisted options
        # and options provided by module
        required_options.difference_update(provided_options)
        required_options.difference_update(blacklisted_options)

        self._config_definition = {'required': required_options,
                                   'default': default_options,
                                   'provided': provided_options,
                                   'blacklisted': blacklisted_options,
                                   #'internal', internal_options,
                                   'excluded': excluded_options,
                                   'iterable': iterable_options
                                   }


    def is_valid(self, modman):

        def custom_cfg_check(options_module):
            if hasattr(options_module, 'check_cfg'):
                return options_module.check_cfg(self)
            else:
                return True

        if not self._config_definition:
            self._load_config_definition(modman)

        cfg_definition      = self._config_definition
        required_options    = cfg_definition['required']
        default_options     = cfg_definition['default']
        provided_options    = cfg_definition['provided']
        blacklisted_options = cfg_definition['blacklisted']

        missing_options = self.missing_options(required_options)
        if missing_options:
            print('Following required options are not present: ')
            for option in missing_options:
                print('  ', repr(option))

            return False

        alien_options = set(self.list_options())

        alien_options.difference_update(required_options)
        alien_options.difference_update([opt[0] for opt in default_options])
        #alien_options.difference_update(internal_options)

        if alien_options:
            print('\nFound following alien options in configuration:')

            for (forbidden, name) in zip([provided_options, blacklisted_options],
                                          ['PROVIDED', 'BLACKLISTED']):

                found_aliens = alien_options.intersection(forbidden)
                if found_aliens:
                    print('\n  Options found in configuration, but %s by '
                          'a centrifuge module:' % name)
                    for option in found_aliens:
                        print('    ', option)

                alien_options.difference_update(found_aliens)

            if alien_options:
                print('\n Options found in configuration, but not specified'
                      'by a module:')
                for option in alien_options:
                    print('  ', option)

            return False

        if not modman.traverse_ancestors(self.get_value('exp_type'),
                                         custom_cfg_check, submodule='options'):
            return False

        return True

        ###alien_options.difference_update([opt[0] for opt in default_options])

##################################################################
#                 ModelParameters class                          #
##################################################################

class ModelParameters:
    """
    Parameters class for the centrifuge simulation.
    """

    def __init__(self, cfg):
        self._cfg            = cfg
        cfg_definition       = cfg._config_definition
        excluded_options     = cfg_definition['excluded']
        self._iterable_options_names = cfg_definition['iterable']

        self._iterable_parameters = {}
        self.iteration            = 0
        self.iterations           = 0

        fns = cfg.get_value('omega2g_fns')
        if fns:
            self._omega2g_fns = \
              {key: MethodType(fn, self) for (key, fn) in fns.items()}

        # read values from cfg and add them to model; if a value is an iterable
        # option, the set_value() method updates self._iterable_parameters
        for (option_name, value) in cfg.iterate_values():
            if not option_name in excluded_options:
                self.set_value(option_name, value)

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

        if hasattr(self, 'duration'):
            if self.duration > 0.0: phase = 'a'
            elif self.fh_duration > 0.0: phase = 'g'
            else: phase = 'd'
            self.set_omega2g_fn(phase)

    def get_iterable_value(self, key):
        if key in self._iterable_parameters:
            return self._iterable_parameters[key]
        else:
             return None

    def set_value(self, key, value):
        """
          Set the value of parameter given as 'key'.
        """

        # Keep self._itarable_parameters up-to-date
        if key in self._iterable_options_names:
            if type(value) in [list, tuple]:
                self._iterable_parameters[key] = value
                value = value[0] # initialize with first value
            else:
                # if previously was a list value and now is an atom, remove it
                if key in self._iterable_parameters:
                    del(self._iterable_parameters[key])

        setattr(self, key, value)

    def set_parameters(self, parameters_dict):
        for (key, value) in parameters_dict.items():
            setattr(self, key, value)

        if key == 'n':
            setattr(self, 'm', 1-1/value)

    def get_parameters(self, parameters):
        params = {}
        for key in parameters:
            params[key] = getattr(self, key)
        return params

    def next_iteration(self):
        """
          Assign the next value of the parameters that were given as type list
        """
        self.iteration += 1

        if self.iteration == self.iterations: return False

        self._set_iteration(self.iteration)

        return True

    def init_iteration(self):
        self.iteration = 0
        self._set_iteration(self.iteration)

    def _set_iteration(self, i):
        """
          Assign i-th value of all parameters supplied (as of type) list
        """
        # values of the i-th iteration are stored at (i-1)-th place in list
        for (key, value) in self._iterable_parameters.items():
            setattr(self, key, value[i])

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

    def find_omega2g(self, t):
        """ Empty definition for omega^2/g function """
        raise ValueError('This is an empty definition and has to be replaced.')

    def set_omega2g_fn(self, fn_type):
        self.find_omega2g = self._omega2g_fns[fn_type]
