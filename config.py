"""
This package implements access to configuration and defines
the ModelParameters class which stores the setting obtained
from the configuration files.
"""

try:
    import ConfigParser as configparser
except:
    import configparser
import numpy as np

##################################################################
#                   Configuration class                          #
##################################################################

def flatten(cfg):
    cfg_dict = cfg._cfg_dict

    def dict_flatten(base_dict, dict_tree):
        for (key, value) in dict_tree.items():
            if type(value) == dict:
                dict_flatten(base_dict, value)
            else:
                base_dict[key] = value

    flattened_cfg = Configuration()
    flattened_cfg._cfg_dict = dict_flatten({}, cfg)

    return flattened_cfg

def parse_value(str_value):
    """
      Given a value as string, tries to convert to it's correspending type.
      May be called recursively in the case of nested structures.
    """
    try:
        raw_value = str_value.strip()

        if "." in raw_value or "e" in raw_value or "E" in raw_value:
            return float(setting)
        elif raw_value[0] == "[" and raw_value[-1] == "]":
            return [parse_value(item) for item in raw_value[1:-1].split(",")]
        elif ((raw_value[0] == "'" and raw_value[-1] == "'")
            or (raw_value[0] == '"' and raw_value[-1] == '"')):
            return raw_value[1:-1]
        elif raw_value = 'True':
            return True
        elif raw_value = 'False':
            return False
        elif raw_value[0] == "(" and raw_value[-1] == ")":
            return tuple([parse_value(item)\
                          for item in raw_value[1:-1].split(",")])
        else:
            return int(raw_value)

    except:
        print('Error:Could not parse setting: ', raw_value, '\nExiting...')
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

    def options_p(self, options, section = None):

        cfg_dict = self._cfg_dict

        if self._preserve_sections_p and (not section
                                          or not section in cfg_dict):
            print('options_p error: Section not present in config: ', section)
            exit(1)

        if preserve_sections_p:
            for option in options:
                if not option in cfg_dict[section]:
                    return False
        else:
             for option in options:
                if not option in cfg_dict:
                    return False
        return True

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

                value = eval_item(raw_value)

                if preserve_sections_p:
                    cfg_dict[section][option]=value
                else:
                    cfg_dict[option]=value

        return self

##################################################################
#                 ModelParameters class                          #
##################################################################

class ModelParameters:
    """
    Parameters of the centrifuge
    """
    def __init__(self, cfg = None):
        if cfg:
            self.register_keys(cfg)
        self.register_key('tspan', np.asarray([], dtype=float))
        self.register_key('omega_fall', np.asarray([], dtype=float))

    def register_key(self, key, value):
        if hasattr(self, key):
            raise Exception("Atrribute '%s' already exists !" % key)
        else:
            self.set_value(key, value)

    def register_keys(self, flattened_cfg):
        for (key, value) in flattened_cfg.items():
            self.register_key(key.lower(), value)

    def set_value(self, key, value):

        setattr(self, key, value)

        if key == 'n':
            if type(value) == list:
                m =  [1-1/n for n in value]
            else:
                m = 1-1/value
            setattr(self, 'm', m)

    def set(self, key, value):
        key_lower = key.lower()

        if not hasattr(self, key_lower):
            raise Exception('Atrribute ''%s'' does not exist !' % key)

        value_type = type(value)
        key_type   = type(getattr(self, key_lower))

        if value_type == key_type:
            setattr(self, key_lower, value)
        elif value_type == int and key_type == float:
            #convert int automatically to float
            setattr(self, key, float(value))
        elif value_type == list:
            for item in value:
                if type(item) == int and key_type == float:
                    pass
                elif not ((type(item) == key_type) or
                          (type(item) == int and key_type == float)):
                    raise ValueError("ModelParameters: key '%s' has wrong type."
                            " Expected type '%s' and got type '%s' of %s"
                            % (key, key_type, value_type, value))
                if value and type(value[0] == int) and (key_type == float):
                    value = [float(item) for item in value]

                setattr(self, key, value)
        else:
            raise ValueError("ModelParameters: key '%s' has wrong type."
                             " Expected type '%s' and got type '%s' of %s"
                             % (key, key_type, value_type, value))

    
