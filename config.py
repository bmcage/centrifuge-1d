"""
This package implements access to configuration and defines
the ModelParameters class which stores the setting obtained
from the configuration files.
"""

try:
    import ConfigParser as configparser
except:
    import configparser
from shared import load_modules, load_experiment_types
from numpy import pi

exptype2modulename = load_experiment_types()

def load_all_options():
    config_options_modules_list = {}

    def traverse(result, module_name):
        if not module_name in result:
            options_module = find_options_module(module_name)

            for parental_module in options_module.PARENTAL_MODULES:
                traverse(result, parental_module)
            result.append(module_name)

    def find_options_modules_names(module_name):
        if module_name in config_options_modules_list:
            return config_options_modules_list[module_name]
        else:
            modules_names_list = []
            traverse(modules_names_list, module_name)
            config_options_modules_list[module_name] = modules_names_list
            return modules_names_list

    return find_options_modules_names

find_options_modules_names = load_all_options()
find_options_module = load_modules('options')

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

    def set_defaults(self):
        cfg_dict = self._cfg_dict

        if not 'exp_type' in self._cfg_dict:
            print('Missing option: configuration does not specify ''exp_type'' '
                  'option. Cannot validate configuration.\nCurrent '
                  'configuration is:')
            self.echo()
            print('\nExiting...')
            exit(1)

        # set options that are missing but have specified default value
        module_name = exptype2modulename(cfg_dict['exp_type'])
        config_modules_list = find_options_modules_names(module_name)
        for config_options_module in config_modules_list:
            current_module = find_options_module(config_options_module)
            options = current_module.CONFIG_OPTIONS['defaults']

            if self._preserve_sections_p:
                for (section, section_content) in options.items():
                    for (option, value) in section_content.items():
                        if not option in cfg_dict:
                            cfg_dict[section][option] = value
            else:
                for (option, value) in options.items():
                    if not option in cfg_dict:
                        cfg_dict[option] = value
        return self

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

    def is_valid(self):

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

        if not 'exp_type' in self._cfg_dict:
            print('Missing option: configuration does not specify ''exp_type'' '
                  'option. Cannot validate configuration.\nCurrent '
                  'configuration is:')
            self.echo()
            print('\nExiting...')
            exit(1)

        module_name = exptype2modulename(self._cfg_dict['exp_type'])
        config_modules_list = find_options_modules_names(module_name)

        alien_options = set(self.list_options())

        for config_module_name in config_modules_list:

            current_module   = find_options_module(config_module_name)
            config_parameters = current_module.CONFIG_OPTIONS

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

        if alien_options:
            print('\nFound following alien options in configuration:')
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
    def __init__(self, cfg, parameters_list):
        self._cfg = cfg
        # rpm->rad.s-1:  omega_radps = (2pi)*omega_rps/60
        self.rpm2radps = lambda x: x * pi/ 30.0
        self._parameters_list = parameters_list
        self._iterable_parameters = []
        self._iteration = 0
        self.first_iteration_p = False

        for param in parameters_list:
            self.set_value(param, cfg.get_value(param))

        self._iterations_count = len(cfg.get_value('duration'))

    def set_value(self, key, value):
        """
          Set the value of parameter given as 'key'.
          Performs the check of the 'value' and conversion to correct units
          when needed. Updates also the dependent variables.
        """
        # Keep self._itarable_parameters up-to-date; if we set a list-type value
        # should be stored, if an atom, should be removed
        if (type(value) == list):
            if not key in self._iterable_parameters:
                self._iterable_parameters.append(key)
        else:
            if key in self._iterable_parameters:
                self._iterable_parameters.remove(key)

        if key in ['omega', 'omega_start', 'omega_end']:
            if type(value) == list:
                setattr(self, key, [self.rpm2radps(omega) for omega in value])
            else:
                setattr(self, key, self.rpm2radps(value))
            return
        elif key in ['omega_fall', 'm']:
            raise ValueError('Parameter %s is not allowed to set directly.'
                             % key)

        setattr(self, key, value)

        if key == 'n':
            if type(value) == list:
                m =  [1-1/n for n in value]
            else:
                m = 1-1/value
            setattr(self, 'm', m)
        elif key == 'r0_fall':
            from math import sqrt

            if type(value) == list:
                setattr(self, key, [sqrt(self.g/r0_fall) for r0_fall in value])
            else:
                setattr(self, key, sqrt(self.g/value))

    def next_iteration(echo):
        """
          Assign the next value of the parameters that were given as type list
        """
        i = self._iteration
        self.first_iteration_p = (i == 0)
        cfg = self._cfg

        for key in self._iterable_parameters:
            setattr(key, cfg.get_value(key)[i])

        self.iteration = i+1

        return (i < self._iterations_count)

    def echo(self, iterable_only=False):
        """
          Print the parameters stored in the model.
          If 'iterable_only' is True, prints only the current values of the
          variables that are iterated by the 'next_iteration()' function.
        """
        if iterable_only:
             parameters = sorted(self._iterable_parameters)
        else:
             parameters = sorted(self._parameters_list)
        print()
        for option in parameters:
            print('%-12s = %s' % (option, getattr(self, option)))

##################################################################
#           Functions on Configuration/ModelParameters           #
##################################################################

def determine_model_options(config_parameters, cfg, cfg_only_options):

    model_options = list(config_parameters['mandatory'])
    for option in cfg_only_options:
        model_options.remove(option)

    model_options.extend(config_parameters['defaults'].keys()
                         + config_parameters['additional'])

    optional_options = config_parameters['optional']
    cfg_missing_options = cfg.missing_options(optional_options)
    for option in optional_options:
        if not option in cfg_missing_options:
            model_options.append(option)

    for (test_fn, dependent_options) in config_parameters['dependent'].values():
        if test_fn(cfg):
            model_options.extend(dependent_options)

    return  model_options
