import numpy as np
from const import INI_DIR, MASKS_DIRNAME

def yn_prompt(question_str):
    while True:
        answ = input(question_str)
        answ = answ.lower()
        if answ in ['', 'y', 'yes']: return True
        if answ in ['n', 'no']: return False

def get_directories(dirs, exp_id, exp_no, tube_no):

    dir_struct = ['exp_base', 'exp_no', 'tube', 'masks']
    dir_values = (exp_id, str(exp_no), 'tube' + str(tube_no), MASKS_DIRNAME)

    def get_dir(dir_type, base_dir):
        k = dir_struct.index(dir_type)
        return base_dir + '/' + '/'.join(dir_values[:k+1]) + '/'

    def resolve_dirs(*dirs):
        results = []
        for dirtype in dirs:
            if dirtype == 'base':
                results.append(INI_DIR)
            elif dirtype == 'exp_base':
                results.append(get_dir('exp_base', INI_DIR))
            elif dirtype == 'exp_no':
                results.append(get_dir('exp_no', INI_DIR))
            elif dirtype in ['data', 'tube']:
                results.append(get_dir('tube', INI_DIR))
            elif dirtype == 'masks':
                results.append(get_dir('masks', INI_DIR))
            elif dirtype == 'search':
                results.append(resolve_dirs('base', 'exp_base', 'exp_no',
                                            'tube'))
            else:
                raise ValueError('Unknown value for get_directories(): '
                                 '{}'.format(dirs))
        return results

    if not dirs: dirs = 'search'

    if type(dirs) == str:
        return resolve_dirs(dirs)[0]
    else:
        return resolve_dirs(*dirs)

def print_by_tube(tube_number, tube_data):
    print('Tube number: ', tube_number)
    for tdata in tube_data:
        print(tdata[0] / 100.)

def make_collector(tubes_numbers):
    by_tube_collection = {tube_no:[] for tube_no in tubes_numbers}
    fifo_collection = []

    def collection(command, data = None, tube_no = None, print_format_fn=print):
        if command == 'collect':
            if tube_no is None:
                fifo_collection.append(data)
            else:
                by_tube_collection[tube_no].append(data)

        elif command in ['print', 'print-by-tube', 'print-fifo']:
            if command == 'print-by-tube':
                for (tube_no, tube_data) in by_tube_collection.items():
                    print('Tube:', tube_no)
                    for value in tube_data:
                        print_format_fn(value)
            else:
                for value in fifo_collection:
                    print_format_fn(value)

        elif command == 'get':
            return (fifo_collection, by_tube_collection)
        else:
            raise ValueError('Collector: Unknown command: %s.\n Valid commans'
                             ' are: "collect", "print", "print-by-tube",'
                             ' "print-fifo" or "get".' % command)
    return collection
