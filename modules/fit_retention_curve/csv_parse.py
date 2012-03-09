
def experiment_identifiers(indexes):
    exp_identifiers = ['exp_type', 'id']
    
    for exp_identifier in exp_identifiers:
        if not exp_identifier in indexes:
            raise("CSV file: missing experiment identifier: '%s'" 
                  % exp_identifier)

    return exp_identifiers

def read_CSV_row(data, exp_type, row, indexes, csv_fields, ini_fields):
    
    experiment_id   = row[indexes['id']]

    if not experiment_id in data:
        data[experiment_id] = {}

    for (descriptor, value) in zip(csv_fields, row):
        if descriptor in ini_fields:
            if not descriptor in data[experiment_id]:
                data[experiment_id][descriptor] = []
            data[experiment_id][descriptor].append(value)

    if not 'exp_type' in data[experiment_id]:
        data[experiment_id]['exp_type'] = row[indexes['exp_type']]
    elif data[experiment_id]['exp_type'] != row[indexes['exp_type']]:
        raise ValueError(experiment_id
                         + ': cannot have different experiment types: ', value, 
                         ' and ' + experiment_data[experiment_id]['exp_type'])

def write2ini(exp_id, exp_id_struct, out_dir):

    exp_type = exp_id_struct['exp_type']
    h        = exp_id_struct['h']
    exp_no   = 1

    for (descriptor, value) in exp_id_struct.items():
        if descriptor in ['exp_type', 'h']: continue

        fout_filename = out_dir + '/experiment_' + str(exp_no) + '.ini'

        fout = open(fout_filename, mode='w', encoding='utf-8')
        fout.write('[experiment-data]\n')
        fout.write("{:8} = '{}'\n".format('exp_type', exp_type))
        fout.write("{:8} = '{}'\n".format('sample_id', descriptor))
        fout.write('{:8} = [{}]\n'.format('h', ', '.join(h)))
        fout.write('{:8} = [{}]\n'.format('u', ', '.join(value)))
        fout.close()

        exp_no = exp_no + 1

    print('Done.')

    
