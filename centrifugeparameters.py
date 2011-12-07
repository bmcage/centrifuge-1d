#
# Copyright (C) 2011-12  Pavol Kišon
# Copyright (C) 2011-12  Benny Malengier
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

#-------------------------------------------------------------------------
#
# Initialization
#
#-------------------------------------------------------------------------

from numpy import arange

class CentrifugeParameters:
    """
    Parameters of the centrifuge                                      
    """
    def get_instance():
        return CentrifugeParameters()

    get_instance = staticmethod(get_instance)

    def __init__(self, preserve_sections_p = False):
        self. preserve_sections_p =  preserve_sections_p

    def register_key(self, section, key, value):
        if self.preserve_sections_p:
            raise Exception("Preserving sections not implemented")
        else:
            if hasattr(self, key):
                raise Exception('Atrribute ''%s''already exists !' % key)
            else:
                setattr(self, key, value)

    def register_keys(self, sections_keys_values_dict):
        for section,keys_values in sections_keys_values_dict.items():
            for key, value in keys_values.items():
                self.register_key(section.lower(), key.lower(), value)

    def set(self, section, key, value):
        key_lower = key.lower()
        if self.preserve_sections_p:
            raise Exception("Preserving sections not implemented")
        else:
            if not hasattr(self, key_lower):
                raise Exception('Atrribute ''%s'' does not exist !' % key)

            if type(value) == type(getattr(self, key_lower)):
                setattr(self, key_lower, value)
            elif type(value) == int and \
                    type(getattr(self, key_lower)) == float:
                #convert int automatically to float
                setattr(self, key, float(value))
            elif type(value) == list:
                key_type = type(getattr(self, key_lower))
                for item in value:
                    if not (type(item) == key_type):
                        print("CentrifugeParameters::WARNING: key has non-string value !"
                              "'%s.%s'" % (section, key))
                        return
                setattr(self, key, value)
            else:
                print("CentrifugeParameters::WARNING: ignoring key with wrong type "
                       "'%s.%s'" % (section, key))


