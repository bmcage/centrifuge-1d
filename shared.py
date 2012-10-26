import numpy as np
from const import INI_DIR, MASKS_DIRNAME, FIGS_DIR

def yn_prompt(question_str):
    while True:
        answ = input(question_str)
        answ = answ.lower()
        if answ in ['', 'y', 'yes']: return True
        if answ in ['n', 'no']: return False

def get_directories(basedir_type, dirs, experiment_info):

    dir_struct = ['exp_base', 'exp_no', 'masks', 'mask']
    dir_values = (experiment_info['exp_id'], str(experiment_info['exp_no']),
                  MASKS_DIRNAME, experiment_info['mask'])

    def get_dir(dir_type, base_dir):
        if dir_type == 'base':
            return base_dir
        elif dir_type in dir_struct:
            k = dir_struct.index(dir_type)
            return base_dir + '/'.join(dir_values[:k+1]) + '/'
        else:
            raise ValueError('Unknown value for get_directories(): '
                                 '{}'.format(dir_type))
    def resolve_dirs(basedir, *dirs):
        results = []
        for dirtype in dirs:
            if dirtype == 'search':
                results.append(resolve_dirs(basedir, 'base', 'exp_base',
                                            'exp_no'))
            elif dirtype == 'data':
                results.append(get_dir('exp_no', basedir))
            else:
                results.append(get_dir(dirtype, basedir))

        return results

    if not dirs: dirs = 'search'

    if basedir_type == 'ini':
        basedir = INI_DIR + '/'
    elif basedir_type == 'figs':
        basedir = FIGS_DIR + '/'
    else:
        raise ValueError('Unrecognized type of basedir_type: ', basedir_type)

    if type(dirs) == str:
        return resolve_dirs(basedir, dirs)[0]
    else:
        return resolve_dirs(basedir, *dirs)

def print_by_tube(tube_number, tube_data):
    print('Tube number: ', tube_number)
    for tdata in tube_data:
        print(tdata[0] / 100.)

##def make_collector(tubes_numbers):
##    by_tube_collection = {tube_no:[] for tube_no in tubes_numbers}
##    fifo_collection = []
##
##    def collection(command, data = None, tube_no = None, print_format_fn=None):
##        if command == 'collect':
##            if tube_no is None:
##                fifo_collection.append(data)
##            else:
##                by_tube_collection[tube_no].append(data)
##
##        elif command in ['print', 'print-by-tube', 'print-fifo']:
##            if command == 'print-by-tube':
##                for (tube_no, tube_data) in by_tube_collection.items():
##                    print('Tube:', tube_no)
##                    for value in tube_data:
##                        print_format_fn(value)
##            else:
##                for value in fifo_collection:
##                    print_format_fn(value)
##
##        elif command == 'get':
##            return (fifo_collection, by_tube_collection)
##        else:
##            raise ValueError('Collector: Unknown command: %s.\n Valid commans'
##                             ' are: "collect", "print", "print-by-tube",'
##                             ' "print-fifo" or "get".' % command)
##    return collection
