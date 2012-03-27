# -*- coding: utf-8 -*-
#
# Centrifuge config
#
# Copyright (C) 2005-2007  Donald N. Allingham
# Copyright (C) 2008-2009  Gary Burton 
# Copyright (C) 2009       Doug Blank <doug.blank@gmail.com>
# Copyright (C) 2009       Benny Malengier <bm@cage.ugent.be>
# Copyright (C) 2011-12    Pavol Kišon     <pk@cage.ugent.be>
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
#


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

#import const

def print_cfg(cfg):
    print()
    for (sec, value) in cfg.items():
        print('section: ', sec)
        print(value)


def flatten_cfg(cfg):
    def flatten(flattened_cfg, cfg):
        for (key, value) in cfg.items():
            if type(value) == dict:
                flatten(flattened_cfg, value)
            else:
                flattened_cfg[key] = value
        return flattened_cfg

    return flatten({}, cfg)

def merge_flattened_cfgs(cfg, *cfgs):
    for acfg in cfgs:
        if acfg:
            cfg.update(acfg)

    return cfg

def merge_cfgs(cfg, *cfgs):
    """
    Merge all following cfgs into 'cfg'; if the same values appear, the last one
    (from last cfg) applies
    """

    def merge_cfg(cfg, ncfg):
        for key in ncfg.keys():
            if key in cfg:
                cfg[key].update(ncfg[key])
            else:
                cfg[key]=dict(ncfg[key])

    if type(cfgs) == list:
        for acfg in cfgs:
            if acfg:
                merge_cfg(cfg, acfg)
    else:
        merge_cfg(cfg, cfgs)

    return cfg

def eval_item(setting):
    """
    Given a value from an ini file, return it in its proper type.
    May be recursively called, in the case of nested structures.
    """
    try:
        setting = setting.strip()

        value = None
        if setting.startswith("'") and setting.endswith("'"):
            value = setting[1:-1]
        elif setting.startswith("[") and setting.endswith("]"):
            list_data = setting[1:-1]
            value = [eval_item(item) for item in list_data.split(",")]
        elif setting == "True":
            value = True
        elif setting == "False":
            value = False
        elif "." in setting or "e" in setting or "E" in setting:
            value = float(setting)
        else:
            value = int(setting)
        return value
    except:
        print('Error:Could not parse setting: ', setting, '\nExiting...')
        exit(1)


def read_cfgs(cfgs_filenames, base_cfg = None, preserve_sections_p=True,
                 cfgs_merge=True):
    """
      Reads config .ini files. If preserve_sections_p is false, removes sections
      and leaves only values from all sections; in case of two sections contain
      the same values, the latter will be used. If cfgs_merge=True, merges all
      the config files into one.
    """
    if not base_cfg:
        base_cfg = {}

    def parser2cfg(parser):
        cfg = base_cfg

        for psection in parser.sections():
            section = psection.lower()

            for option in parser.options(psection):
                setting = parser.get(psection, option).strip()

                value = eval_item(setting)

                if preserve_sections_p:
                    if section in cfg:
                        cfg[section][option]=value
                    else:
                        cfg[section]={option: value}
                else:
                    cfg[option]=value

        return cfg

    if not (type(cfgs_filenames) == list):
        cfgs_filenames = [cfgs_filenames]

    if cfgs_merge:
        parser = configparser.ConfigParser()

        read_files = parser.read(cfgs_filenames)

        parsers = [parser]

    else:
        read_files = []
        parsers = []

        for filename in cfgs_filenames:
            parser = configparser.ConfigParser()
            parsers.append(parser)

            if parser.read(filename):
                read_files.append(filename)

    if (len(read_files) != len(cfgs_filenames)):
        print('From expected files: ', str(cfgs_filenames),
             'were successfully parsed only: ', str(read_files))

    cfgs = [parser2cfg(parser) for parser in parsers]

    return cfgs

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
            #print(cfg)
        self.register_key('tspan', np.asarray([], dtype=float))
        self.register_key('omega_fall', np.asarray([], dtype=float))

    def register_key(self, key, value):
        if hasattr(self, key):
            raise Exception("Atrribute '%s' already exists !" % key)
        else:
            setattr(self, key, value)

    def register_keys(self, flattened_cfg):
        for (key, value) in flattened_cfg.items():
            self.register_key(key.lower(), value)

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
            raise ValueError("ModelParameters: key '%s' has wrong type. Expected"
                             " type '%s' and got type '%s' of %s"
                             % (key, key_type, value_type, value))

    
