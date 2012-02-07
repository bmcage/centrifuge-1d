from modules.direct_saturated.csv_parse import filtered_fields

def skip(exp_type, row, indexes):
    if exp_type == 'ish':
        # only those rows where the difference between the
        # starting and ending soil heigh has less than 5% difference
        sh0 = float(row[indexes['sh0']])
        sh1 = float(row[indexes['sh1']])
        return (abs((sh0 - sh1) / sh0) > 0.05)
    else:
        return False

def CSV2ini_fields(field_name, field_value)
    return (field_name, field_value)
