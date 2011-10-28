#
# Copyright (C) 2011-12  Pavol Kišon
# Copyright (C) 2009-12  Benny Malengier
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

from sys import argv as sysargv
from config import DEFAULT_PARAMETERS, ConfigManager
from numpy import arange

from centrifugeparameters import CentrifugeParameters
from auxiliaryfunctions   import parse_centrifuge_input

[inifiles, outputdir, saveconfigname] = \
    parse_centrifuge_input(sysargv)

def main():
    cfgmngr = ConfigManager.get_instance()
    model   = cfgmngr.read_configs(merge_output_p = True, \
                                   preserve_sections_p = False, \
                                   filenames = inifiles, \
                                   defaults = [DEFAULT_PARAMETERS], \
                                   saveconfigname)
    model.register_key('experiment', 'tspan', \
                   arange(model.t_start, model.t_end, model.t_step))

if __name__ == "__main__":
    main()
