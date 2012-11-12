# -*- coding: utf-8 -*-
from os import sep, getcwd

#-------------------------------------------------------------------------
# Paths
#-------------------------------------------------------------------------
PROGRAM_DIR = getcwd()
DATA_DIR    = PROGRAM_DIR + sep + 'data'
CSV_DIR     = DATA_DIR    + sep + 'csv'
INI_DIR     = DATA_DIR    + sep + 'datafiles'
FIGS_DIR    = DATA_DIR    + sep + 'datafiles'

#-------------------------------------------------------------------------
# Default file and directory names
#-------------------------------------------------------------------------
DEFAULTS_ININAME = 'defaults.ini'
CONSTANTS_ININAME = 'constants.ini'
PLOTSTYLE_ININAME = 'plotstyles.ini'
MASKS_DIRNAME    = 'masks'
# MEASUREMENTS_NAMES are of type: name:option_name
MEASUREMENTS_NAMES = {'MI': 'wl1', 'MO': 'wl_out', 'GC': 'gc1', 'RM': 'rm1',
                      'theta': 'theta'}
#-------------------------------------------------------------------------
#-------------------------------------------------------------------------
