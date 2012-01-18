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
                    if type(item) == int and key_type == float:
                        pass
                    elif not ((type(item) == key_type) or
                              (type(item) == int and key_type == float)):
                        print("CentrifugeParameters::WARNING: key '%s.%s' has non-expected type: %s, expected type: %s" 
                              % (section, key, item, type(item), key_type))
                        return

                if value and type(value[0] == int) and (key_type == float):
                    value = [float(item) for item in value]

                setattr(self, key, value)
            else:
                print("CentrifugeParameters::WARNING: ignoring key with wrong type "
                       "'%s.%s'" % (section, key))


