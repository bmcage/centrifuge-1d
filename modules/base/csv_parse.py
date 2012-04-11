def skip(exp_type, row, indexes):
    """
      Check every 'row' of a .csv file. If the row is not wanted, set return
      flag to 'False' and row will be ommited from the resulting .ini file.
      Otherwise set return flag to 'True'.
    """
    return False

def CSV2ini_fields(field_name, field_value):
    return (field_name, field_value)
