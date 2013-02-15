def skip(row):
    """
      Check every 'row' of a .csv file. If the 'row' should be discarded, set
      return flag to 'False' (so it will be ommited from the resulting .ini
      file). Otherwise set return flag to 'True'. 'Row' is of type 'namedtuple'.
    """
    return False

def CSV2ini_fields(field_name, field_value):
    return (field_name, field_value)
