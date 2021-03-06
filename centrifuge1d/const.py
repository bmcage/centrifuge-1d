# -*- coding: utf-8 -*-
from __future__ import print_function
from os import sep
from os.path import dirname, normpath
from pickle import HIGHEST_PROTOCOL

#-------------------------------------------------------------------------
# Paths
#-------------------------------------------------------------------------
PROGRAM_DIR = dirname(__file__)  # location of const.py !
PROTO_DIR   = PROGRAM_DIR + sep + 'proto' + sep
DATA_DIR    = normpath(PROGRAM_DIR + sep + '..' + sep + 'data')
print('Program and data directory set to:', PROGRAM_DIR, DATA_DIR)
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
DIFFPROG = "diff"
