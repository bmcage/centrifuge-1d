def skip(exp_type, row, indexes):
    return False

def CSV2ini_fields(field_name, field_value):
    return (field_name, field_value)

def experiment_identifiers(indexes):
    exp_identifiers = ['exp_type', 'id', 'exp_no', 'tube_no']
    for exp_identifier in exp_identifiers:
        if not exp_identifier in indexes:
            raise("CSV file: missing experiment identifier: '%s'" 
                  % exp_identifier)

    return exp_identifiers

def read_CSV_row(data, exp_type, row, indexes, csv_fields, ini_fields):

    experiment_id   = row[indexes['id']]
    experiment_no   = row[indexes['exp_no']]
    tube_number     = row[indexes['tube_no']]

    if not experiment_id in data:
        data[experiment_id]={}
    if not experiment_no in data[experiment_id]:
        data[experiment_id][experiment_no]={}

    experiment_data = data[experiment_id][experiment_no]
    if not tube_number in experiment_data:
        experiment_data[tube_number]  = {k: [] for k in ini_fields}

    for (descriptor, value) in zip(csv_fields, row):
        if descriptor in ini_fields:
            experiment_data[tube_number][descriptor].append(value)

    if not 'exp_type' in experiment_data[tube_number]:
        experiment_data[tube_number]['exp_type'] = row[indexes['exp_type']]
    elif experiment_data[tube_number]['exp_type'] != row[indexes['exp_type']]:
        raise ValueError(experiment_id +': experiment ' + experiment_no 
                         + ' - cannot have different experiment'
                         + ' types: ', value, ' and '
                         + experiment_data[tube_number]['exp_type'])

def write2ini(exp_id, exp_id_struct, out_dir):
    for (experiment_no, exp_no_struct) in exp_id_struct.items():
        print('        Processing experiment number: %s... ' % experiment_no, end="")

        for (tube_no, base_data) in exp_no_struct.items():

                fout_filename = ''.join([out_dir, '/experiment_', experiment_no,
                                         '-filter', tube_no, '.ini'])

                fout = open(fout_filename, mode='w', encoding='utf-8')
                fout.write('[experiment-data]\n')

                for (descriptor, value) in base_data.items():
                    if descriptor == 'exp_type':
                        fout.write("{:8} = '{}'\n".format(descriptor, value))
                    else:
                        fout.write('{:8} = [{}]\n'.format(descriptor,
                                                          ', '.join(value)))
                fout.close()
        print('Done.')
