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

#import os
#import time
try:
    import ConfigParser as configparser
except:
    import configparser
  #import errno
#from gettext import gettext as _

#import const

def base_cfg():
    base = { 
            'general': {'g': 981., 'debugging': False},
    'starting-filter': {'d1': 0., 'ks1': -1.0 },
               'soil': {'n': 2.81, 'gamma': 0.0189, 'ks': 2.4e-5,
                        'l': 10.0, 'porosity': 0.4, 'v': 1.0},
      'ending-filter': {'d2': 0., 'ks2': -1. },
              'fluid': {'viscosity': 1.0, 'density': 1.0,
                        's1_0': 0.1, 's2_0': 0.2, 'pc0': 1.0e5 },
         'centrifuge': {'r0': 30.0, 'l0_in': 2.0, 'l0_out': 4.0,
                        'd': 4.0, 'deceleration_duration': 0.},
         'experiment': {'exp_type': '',
                        't_start': 0.0, 't_end': 2000.0, 't_step': 200.0,
                        'omega_start': 0.0, 'omega': 35.0, 'omega_gamma': 0.5,
                        'omega_end': 0.0,
                        'inverse_data_filename': '', 'data_type': 0},
     'discretization': {'inner_points': 80, 'first_point_offset': 80.0, 'dtype': 1,
                        'percent_in_saturation': 40.0,
                        'approximation_type': 5, 'mb_epsilon': 1e-5}
    }
    return base

def print_cfg(cfg):
    print()
    for (sec, value) in cfg.items():
        print('section: ', sec)
        print(value)

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
            merge_cfg(cfg, acfg)
    else:
        merge_cfg(cfg, cfgs)

    return cfg

DEFAULT_PARAMETERS = base_cfg()

inverse_cfg = {'inverse_data': {'duration': -1.0, 'h0': -1.0, 'h1': -1.0,
                                'r0': -1.0, 'length': -1.0, 'omega': -1.0,
                                'porosity': -1.0,
                                'd1': 0.0, 'ks1': -1.0, 'd2': 0.0, 'ks2': -1.0,
                                'exp_type': ''}}

DEFAULT_DATA_PARAMETERS = merge_cfgs(base_cfg(), inverse_cfg)

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


# def parser2register(parser, register_fns):
#         """ 
#         Loads data from parser into model.
#         """

#         for psection in parser.sections():
#             section = sec.lower()

#             for option in parser.options(psection):
#                 setting = parser.get(psection, option).strip()
#                 value = eval_item(setting)
#                 register_fn(section, option, value)


def read_configs(cfgs_filenames, preserve_sections_p=True, cfgs_merge=True):
    """
      Reads config .ini files. If preserve_sections_p is false, removes sections
      and leaves only values from all sections; in case of two sections contain
      the same values, the latter will be used. If cfgs_merge=True, merges all
      the config files into one.
    """

    def parser2cfg(parser):
        cfg = {}

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


    # def read_configs_old(self, merge_output_p = True, preserve_sections_p = False, \
        #                  filenames = [const.INIFILE_DEFAULT], \
        #                  defaults = [DEFAULT_PARAMETERS], \
        #                  saveconfigname = ''):
        # """
        # Reads .ini files and returns them as an object(s)
        # """

        # if merge_output_p:
        #     model = centrifugeparameters.CentrifugeParameters(preserve_sections_p)
        #     for default in defaults:
        #         model.register_keys(default)

        #     parser = configparser.ConfigParser()
        #     if saveconfigname:
        #         # save the first filename that is expected to be config
        #         # (so save it without comments)
        #         raise Exception('Saving config file is not implemented !')

        #     for filename in filenames:
        #          if  not (filename and os.path.exists(filename)):
        #              raise ValueError('The file %s cannot be found' % filename)
        
        #          parser.read(filename)
        #     self._parser2model(parser, model)
        #     return model
        # else:
        #     raise Exception("Multiple outputs are not implemented")

    
