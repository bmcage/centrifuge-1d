# -*- coding: utf-8 -*-
from os import sep, getcwd
#from platform import python_version_tuple

#-------------------------------------------------------------------------
#
# Paths
#
#-------------------------------------------------------------------------
PROGRAM_DIR = getcwd()
DATA_DIR    = PROGRAM_DIR + sep + 'data'
CSV_DIR     = DATA_DIR    + sep + 'csv'
INI_DIR     = DATA_DIR    + sep + 'datafiles'
FIGS_DIR    = DATA_DIR    + sep + 'datafiles'

# Default names
DEFAULTS_ININAME = 'defaults.ini'
CONSTANTS_ININAME = 'constants.ini'
PLOTSTYLE_ININAME = 'plotstyles.ini'
MASKS_DIRNAME    = 'masks'
# MEASUREMENTS_NAMES are of type: name:option_name
MEASUREMENTS_NAMES = {'MI': 'wl1', 'MO': 'wl_out', 'GC': 'gc1', 'RM': 'rm1',
                      'theta': 'theta'}
#-------------------------------------------------------------------------
#
# About box information
#
#-------------------------------------------------------------------------
# URL_HOMEPAGE    = "http://cage.ugent.be/~bm/progs"

# PROGRAM_NAME   = "Centrifuge"
# VERSION        = "0.1.1a"
# #COPYRIGHT_MSG  = u"\u00A9 2009 Pavol Kišon" \
# #                  u"\u00A9 2009 Benny Malengier"
# COMMENTS       = "Centrifuge is a simulator of an water infiltration "\
#                  "into soil sample put in a centrifuge"
# PYTHONVERSION  = python_version_tuple()
