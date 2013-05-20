"""
This package implements access to configuration and defines
the ModelParameters class which stores the setting obtained
from the configuration files.
"""

from __future__ import print_function, division

import numpy as np

try:
    import ConfigParser as configparser
except:
    import configparser
from const import DEFAULTS_ININAME, CONSTANTS_ININAME, MEASUREMENTS_ININAME
from shared import parse_value, get_directories, yn_prompt
from os import listdir, path, makedirs
import sys
from sys import modules as sysmodules
from types import MethodType
from modules.shared.measurements import Measurements

##################################################################
#                   ModulesManager class                         #
##################################################################

def get_ancestors(options_module):
    return options_module.PARENTAL_MODULES

class ModulesManager():
    def __init__(self):
        available_modules = listdir('modules')
        for name in ('__init__.py', '__init__.pyc', '__pycache__'):
            if name in available_modules: available_modules.remove(name)

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
                print('INFO: Information about module "%s" could not be '
                      'loaded. Module loading is skipped.'  % module_name)

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
            module = find_module(module_name, submodule='info', not_found=None)
            if module is None: continue
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

    def find_module(self, modname_or_exptype, submodule='', not_found=None):
        """
          Return an module determined by it's name or the type of experiment.
          Loading of modules is done only on demand.
        """

        if not modname_or_exptype in self._modules_names:
            print('\nFind module error: Unknown module name or experiment '
                  'type: "'+ modname_or_exptype + '".\tSkipping...')
            return not_found

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

        if key in cfg_dict:
            del cfg_dict[key]

    def iterate_values(self, section = None):
        if self._preserve_sections_p and section:
            print('cfg:iterate_values: preserving sections is not implemented.')
            exit(1)

        return self._cfg_dict.items()

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

    def _group(self):
        measurements = Measurements()
        measurements.read(self)

        self.set_value('measurements', measurements)

    def adjust_cfg(self, modman):

        def module_preadjust_cfg(module):
            if hasattr(module, 'prior_adjust_cfg'):
                module.prior_adjust_cfg(self)
            return True

        def module_adjust_cfg(module):
            if hasattr(module, 'adjust_cfg'):
                module.adjust_cfg(self)
            return True

        # Before processing custom adjustment set provided options
        # (by default to None)
        cfg_dict = self._cfg_dict
        for provided_option in self._config_definition['provided']:
            cfg_dict[provided_option] = None

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
            if sys.version_info[0] < 3:
                try:
                    parser.read(fname)
                except Exception as E:
                    print(E)
                    exit(0)
            else:
                try:
                    parser.read(fname)
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

    def load_definition(self, modman):

        if not self._preserve_sections_p:
            section = 'foo'
        else:
            raise NotImplementedError('Section preserving is not implemented.')

        cfg_dict = self._cfg_dict

        parsed_options = {'required': set([]),
                          'default': [],
                          'provided': set([]),
                          'blacklisted': set([]),
                          'internal': set([]),
                          'excluded': set([]),
                          'iterable': set([])}

        name2id = {'PROVIDE_OPTIONS': 'provided',
                   'BLACKLIST_OPTIONS': 'blacklisted',
                   'INTERNAL_OPTIONS': 'internal',
                   'EXCLUDE_FROM_MODEL': 'excluded',
                   'OPTIONS_ITERABLE_LISTS': 'iterable'}

        def parse_options(options):
            required_options  = []
            default_options   = []

            while options:
                defaults = []
                dependent_options = []

                # options classification: required, defaults, dependent
                for option in options:
                    option_type = type(option)
                    if option_type == str:
                        required_options.append(option)
                    elif (option_type == list) or (option_type == tuple):
                        if option == []: continue

                        if hasattr(option[0], '__call__'):
                            dependent_options.append(option)
                        else:
                            defaults.append(option)
                    elif hasattr(option, '__call__'):
                        dependent_options.append(option)
                    else:
                        print('Unknown option type:', option)
                        exit(1)

                # setting the defaults values if option not present
                default_options.extend(defaults)
                for (option, value) in defaults:
                    if not option in cfg_dict:
                        self.set_value(option, value, section)

                # all supplied  options were processed
                # build new options group by parsing the dependent options
                options = []

                # process dependent options
                for option in dependent_options:
                    if hasattr(option, '__call__'):
                        new_opts = option(self)
                    elif option[0](self):
                        new_opts = option[1]
                    elif len(option) > 2:
                        new_opts = option[2]
                    else:
                        continue

                    if type(new_opts) in (list, tuple):
                        options.extend(new_opts)
                    else:
                        options.append(new_opts)

            return (required_options, default_options)

        def read_options_definition(options_module):
            # options_name = 'INTERNAL_OPTIONS'
            for (options_name) in ('CONFIG_OPTIONS', 'PROVIDE_OPTIONS',
                                   'BLACKLIST_OPTIONS', 'EXCLUDE_FROM_MODEL',
                                   'OPTIONS_ITERABLE_LISTS'):

                if not hasattr(options_module, options_name): continue

                (required, defaults) = \
                  parse_options(getattr(options_module, options_name))

                if options_name == 'CONFIG_OPTIONS':
                    parsed_options['required'].update(required)
                    parsed_options['default'].extend(defaults)
                else:
                    if defaults:
                        print("Options in '" + options_name + "' cannot have a "
                              "default value:\n",
                              getattr(options_module, options_name),
                              '\nExiting...')
                        exit(1)

                    parsed_options[name2id[options_name]].update(required)

                    # from required options remove blacklisted options
                    # and options provided by module
                    if options_name in ['PROVIDE_OPTIONS', 'BLACKLIST_OPTIONS']:
                        parsed_options['required'].difference_update(required)

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

        self._config_definition = parsed_options

    def is_valid(self, modman, verbose=True):

        def custom_cfg_check(options_module):
            if hasattr(options_module, 'check_cfg'):
                return options_module.check_cfg(self)
            else:
                return True

        if not self._config_definition:
            self.load_definition(modman)

        cfg_definition      = self._config_definition
        required_options    = cfg_definition['required']
        default_options     = cfg_definition['default']
        provided_options    = cfg_definition['provided']
        blacklisted_options = cfg_definition['blacklisted']

        missing_options = self.missing_options(required_options)
        if missing_options:
            if verbose:
                print('\nThe full configuration is:')
                self.echo()

            print('\nConfiguration is NOT VALID.'
                  '\nFollowing required options are not present: ')
            for option in missing_options:
                print('  ', repr(option))
            print()

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
                print('\n  Options found in configuration, but not specified '
                      'by a module:')
                for option in alien_options:
                    print('   ', option)

            if not yn_prompt('\nAlien options found. Do you wish to '
                             'continue? [Y/n]: ', default='y'):
                return False

        if not modman.traverse_ancestors(self.get_value('exp_type'),
                                         custom_cfg_check, submodule='options'):
            return False

        return True

    def set_parameters(self, parameters_dict, section=None):
        if section: raise ValueError('Sections not supported')

        self._cfg_dict.update(parameters_dict)

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

        if hasattr(self, 'SC'):
            # set transformation fns in case we are unsaturated and use
            #  parameterstrasformation
            if hasattr(self, 'transform_params') and self.transform_params:
                max_value = 1e150 # hack - see base_inverse.options.adjust_cfg()
                self.SC.add_transformations_fns(self._transform,
                                                self._untransform, max_value)

    def _get_iterable_value(self, key, not_found=None):
        if key in self._iterable_parameters:
            return self._iterable_parameters[key]
        else:
            return not_found

    def get_value(self, key, not_found=None):
        if hasattr(self, key):
            return getattr(self, key)
        else:
            return not_found

    def set_value(self, key, value):
        """
          Set the value of parameter given as 'key'.
        """

        # Keep self._itarable_parameters up-to-date
        if key in self._iterable_options_names:
            if type(value) in [list, tuple, np.ndarray]:
                self._iterable_parameters[key] = value
                value = value[0] # initialize with first value
            else:
                # if previously was a list value and now is an atom, remove it
                if key in self._iterable_parameters:
                    del(self._iterable_parameters[key])

        setattr(self, key, value)

    def set_parameters(self, parameters_dict):
        if hasattr(self, 'SC'):
            self.SC.set_parameters(parameters_dict)

        for (key, value) in parameters_dict.items():
            if hasattr(self, key): setattr(self, key, value)

    def get_parameters(self, parameters):
        return {key: getattr(self, key)
                for key in parameters if hasattr(self, key) }

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

    def get_measurements_nr(self):
        meas_nr = 0
        for meas_name in ('measurements_times', 'measurements_dec_times',
                          'measurements_fh_times'):
            value = self._get_iterable_value(meas_name)
            if value:
                for meas in value: meas_nr += len(meas)

        return meas_nr

##################################################################
#                   Load model function                          #
##################################################################

def load_configuration(experiment_info, include_global_constants=True):
    (search_dirs, data_dir, masks_dir) = \
      get_directories('ini', ['search', 'data', 'masks'], experiment_info)

    defaults_files = \
      [fname for fname in [sdir + DEFAULTS_ININAME for sdir in search_dirs]
       if path.exists(fname)]

    measurements_files = [data_dir + MEASUREMENTS_ININAME]

    # Read experiment default.ini files and measurements.ini
    # Handle also case when 'csv_file' option is set
    cfg_files = defaults_files + measurements_files
    cfg = Configuration().read_from_files(*cfg_files)

    csv_datafilename = cfg.get_value('csv_file', not_found=None)
    if csv_datafilename:
        cfg.read_from_files(data_dir + csv_datafilename + '.ini')

    # Read masks files
    masks_files = []
    masks = experiment_info['mask']

    if masks:
        for mask_name in masks:
            mask_filename = masks_dir + mask_name + '.ini'
            if not path.exists(mask_filename):
                print('Mask file "{}" does not exist in expected location:'
                      '\n{}.'.format(mask_name, masks_dir))
                if not yn_prompt('Do you wish to continue without applying '
                                 'the mask? [Y/n]: '):
                    exit(0)

            masks_files.append(mask_filename)

        cfg.read_from_files(*masks_files)

    # Handle CONSTANTS inifiles
    if include_global_constants:
        constants_files = \
          [fname for fname in [sdir + CONSTANTS_ININAME for sdir in search_dirs]
           if path.exists(fname)]

        if constants_files:
            consts_cfg = Configuration().read_from_files(*constants_files)

            tube_no = cfg.get_value('tube_no')
            cfg.set_parameters(consts_cfg.get_value('tubes')[tube_no])

    return cfg

def load_model(experiment_info, display_only=False, validate=True,
               modman = None):

    # Process also global constants
    cfg = load_configuration(experiment_info, include_global_constants=True)

    if display_only:
        header = ("Configuration file of experiment '{}' number {:d}"
                  .format(experiment_info['exp_id'],
                          experiment_info['exp_no']))
        print("\n", header, '\n', len(header) * '-')
        cfg.echo()

        return None

    if modman is None:
        modman = ModulesManager()

    if validate and (not cfg.is_valid(modman, verbose=True)): exit(1)

    cfg.adjust_cfg(modman)

    model = ModelParameters(cfg)

    model.experiment_info = experiment_info

    return model
