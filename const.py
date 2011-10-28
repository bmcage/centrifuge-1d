#
# Copyright (C) 2009-12  Benny Malengier
# Copyright (C) 2011-12  Pavol Kišon
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

INIFILE_DEFAULT = 'centrifuge-default.ini'
INIFILE_DEFAULT_EXT = 'centrifuge-default-extended.ini'

#-------------------------------------------------------------------------
#
# Options Constants
#
#-------------------------------------------------------------------------

LONGOPTS = ['inifiles', 'outputdir', 'saveconfig']
SHORTOPTS = 'i:o:s' 
