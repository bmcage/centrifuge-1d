"""
Provides constants for other modules

To access:
import const
const.VALUE
"""

#-------------------------------------------------------------------------
#
# Standard python modules
#
#-------------------------------------------------------------------------
import os
import platform

#-------------------------------------------------------------------------
#
# Website
#
#-------------------------------------------------------------------------
URL_HOMEPAGE    = "http://cage.ugent.be/~bm/progs"

#-------------------------------------------------------------------------
#
# paths
#
#-------------------------------------------------------------------------
USER_HOME = os.path.expanduser('~') 
PROGRAM_DIR = os.getcwd()
DATA_DIR  = PROGRAM_DIR + os.sep + 'centrifuge'

# dirs that need to be created
USER_DIRLIST = (DATA_DIR,)

#-------------------------------------------------------------------------
#
# About box information
#
#-------------------------------------------------------------------------
PROGRAM_NAME   = "Centrifuge"
VERSION        = "0.0.1a"
#COPYRIGHT_MSG  = u"\u00A9 2009 Pavol Kišon" \
#                  u"\u00A9 2009 Benny Malengier"
COMMENTS       = "Centrifuge is a simulator of an water infiltration "\
                 "into soil sample put in a centrifuge"
PYTHONVERSION  = platform.python_version_tuple()

#-------------------------------------------------------------------------
#
# Constants
#
#-------------------------------------------------------------------------

INIFILE_DEFAULT = 'inifiles-defaults/parameters_default.ini'
INIFILE_DEFAULT_EXT = 'inifiles-defaults/parameters_default_extended.ini'

#-------------------------------------------------------------------------
#
# Options Constants
#
#-------------------------------------------------------------------------

LONGOPTS = ['inifiles', 'outputdir', 'saveconfig']
SHORTOPTS = 'i:o:s' 
