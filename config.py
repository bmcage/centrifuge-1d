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

#---------------------------------------------------------------
#
# System imports
#
#---------------------------------------------------------------
import os
import time
import configparser
import errno
from gettext import gettext as _

#---------------------------------------------------------------
#
# Centrifuge imports
#
#---------------------------------------------------------------
import const, centrifugeparameters

#---------------------------------------------------------------
#
# Constants
#
#---------------------------------------------------------------
DEFAULT_PARAMETERS = \
    {\
            'general': {'g': 981., 'debugging': False}, \
               'soil': {'n': 2.81, 'gamma': 0.0189, 'Ks': 2.4e-5, \
                        'porosity': 0.4, 'V': 1.0}, \
              'fluid': {'viscosity': 1.0, 'density': 1.0, \
                        's1_0': 0.1, 's2_0': 0.2, 'pc0': 1.0e5 }, \
         'centrifuge': {'r0': 30.0, 'L': 10.0, 'l0_in': 2.0, 'l0_out': 4.0,
                        'd': 4.0}, \
         'experiment': {'t_start': 0.0, 't_end': 2000.0, 't_step': 200.0, \
                        'omega_start': 0.0, 'omega': 35.0, 'omega_gamma': 0.5,
                        'inverse_data_filename': '', 'inverse_data_type': 0}, \
     'discretization': {'inner_points': 80, 'first_point_offset': 80.0, 'dtype': 1, \
                        'percent_in_saturation': 40.0, \
                        'approximation_type': 5, 'mb_epsilon': 1e-5}}

DEFAULT_DATA_PARAMETERS = \
    {'inverse_data': {'duration': -1.0, 'h0': -1.0, 'h1': -1.0,
                      'length': -1.0, 'omega': -1.0, 'porosity': -1.0}}

#---------------------------------------------------------------
#
# Local functions
#
#---------------------------------------------------------------
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

#---------------------------------------------------------------
#
# Classes
#
#---------------------------------------------------------------
class ConfigManager(object):
    """
    Class to construct the singleton CONFIGMAN where all 
    settings are stored.
    """
    __instance = None
    
    def get_instance():
        """ Use this function to get the instance of the ConfigManager 
        that will work on inifile
        """
        if ConfigManager.__instance is None:
            ConfigManager.__instance = 1 # Set to 1 for __init__()
            ConfigManager.__instance = ConfigManager()
        return ConfigManager.__instance
    get_instance = staticmethod(get_instance)
    
    def __init__(self):
        """ 
        Configure manager constructor loads filenames into 'model' with default
        values 'defaults'
        """
        if ConfigManager.__instance is not 1:
            raise Exception("This class is a singleton. ", \
                            "Use the get_instance() method")

    def read_configs(self, merge_output_p = True, preserve_sections_p = False, \
                         filenames = [const.INIFILE_DEFAULT], \
                         defaults = [DEFAULT_PARAMETERS], \
                         saveconfigname = ''):
        """
        Reads .ini files and returns them as an object(s)
        """

        if merge_output_p:
            model = centrifugeparameters.CentrifugeParameters(preserve_sections_p)
            for default in defaults:
                model.register_keys(default)

            parser = configparser.ConfigParser()
            if saveconfigname:
                # save the first filename that is expected to be config
                # (so save it without comments)
                raise Exception('Saving config file is not implemented !')

            for filename in filenames:
                 if  not (filename and os.path.exists(filename)):
                     raise ValueError('The file %s cannot be found' % filename)
        
                 parser.read(filename)
            self._parser2model(parser, model)
            return model
        else:
            raise Exception("Multiple outputs are not implemented")

    def _parser2model(self, parser, model):
        """ 
        Loads data from parser into model.
        """

        for sec in parser.sections():
            section = sec.lower()
            
            for option in parser.options(sec):
                setting = parser.get(sec, option).strip()
                value = eval_item(setting)
                model.set(section, option, value)
