# -*- coding: utf-8 -*-
from os import sep, getcwd
from pickle import HIGHEST_PROTOCOL

#-------------------------------------------------------------------------
# Paths
#-------------------------------------------------------------------------
PROGRAM_DIR = getcwd() # Program directory
DATA_DIR    = PROGRAM_DIR + sep + 'data'
CSV_DIR     = DATA_DIR    + sep + 'csv'
INI_DIR     = DATA_DIR    + sep + 'datafiles'
FIGS_DIR    = DATA_DIR    + sep + 'datafiles'

#-------------------------------------------------------------------------
# Default file and directory names
#-------------------------------------------------------------------------
DEFAULTS_ININAME   = 'defaults.ini'
CONSTANTS_ININAME  = 'constants.ini'
MEASUREMENTS_ININAME  = 'measurements.ini'
PLOTSTYLE_ININAME  = 'plotstyles.ini'
MASKS_DIRNAME      = 'masks'
DUMP_DATA_FILENAME = 'data.dat'

#-------------------------------------------------------------------------
# Plotting
#-------------------------------------------------------------------------
# Dumping data format: This specifies the format used for data dumping.
# For compatibility with python 2.x the value of '2' is used. In general,
# if compatibility is not needed, the HIGHEST_PROTOCOL value could be used
DUMP_DATA_VERSION=2
