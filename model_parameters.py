from numpy import arange

class ModelParameters:
    """
    Parameters of the centrifuge
    """
    def __init__(self, cfg = None):
        if cfg:
            self.register_keys(cfg)

    def register_key(self, key, value):
        if hasattr(self, key):
            raise Exception("Atrribute '%s' already exists !" % key)
        else:
            setattr(self, key, value)

    def register_keys(self, cfg):
        for (section, keys_values) in cfg.items():
            for (key, value) in keys_values.items():
                self.register_key(key.lower(), value)

    def set(self, key, value):
        key_lower = key.lower()

        if not hasattr(self, key_lower):
            raise Exception('Atrribute ''%s'' does not exist !' % key)

        value_type = type(value)
        key_type   = type(getattr(self, key_lower))

        if value_type == key_type:
            setattr(self, key_lower, value)
        elif value_type == int and key_type == float:
            #convert int automatically to float
            setattr(self, key, float(value))
        elif value_type == list:
            for item in value:
                if type(item) == int and key_type == float:
                    pass
                elif not ((type(item) == key_type) or
                          (type(item) == int and key_type == float)):
                     raise ValueError("ModelParameters: key '%s' has wrong type."
                             " Expected type '%s' and got type '%s' of %s"
                             % (key, key_type, value_type, value))
                if value and type(value[0] == int) and (key_type == float):
                    value = [float(item) for item in value]

                setattr(self, key, value)
        else:
            raise ValueError("ModelParameters: key '%s' has wrong type. Expected"
                             " type '%s' and got type '%s' of %s"
                             % (key, key_type, value_type, value))


