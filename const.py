"""
Provides constants for other modules

To access:
import const
const.VALUE
"""
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
INI_DIR     = DATA_DIR    + sep + 'inifiles'
FIGS_DIR    = DATA_DIR    + sep + 'figures'

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
