def skip(exp_type, row, indexes):
    if exp_type in ['di', 'dih']:
        # only those rows where the difference between the
        # starting and ending soil heigh has less than 5% difference
        # and the amount of outspelled water is non-negative
        l0 = float(row[indexes['l0']])
        l1 = float(row[indexes['l1']])
        return abs((l0 - l1) / l0) > 0.05
    else:
        return False

def CSV2ini_fields(field_name, field_value):
    return (field_name, field_value)
