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
This package implements access to configuration.
"""

try:
    import ConfigParser as configparser
except:
    import configparser

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

def merge_cfgs(cfg, cfgs):
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
                    cfg[setting]=value

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

    
