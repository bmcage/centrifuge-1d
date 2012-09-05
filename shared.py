import numpy as np
from const import INI_DIR, MASKS_DIRNAME, DEFAULTS_ININAME, CONSTANTS_ININAME

def yn_prompt(question_str):
    while True:
        answ = input(question_str)
        answ = answ.lower()
        if answ in ['', 'y', 'yes']: return True
        if answ in ['n', 'no']: return False

def get_default_ini_filename(initype):
    if initype == 'default':
        return DEFAULTS_ININAME
    elif initype == 'constants':
        return CONSTANTS_ININAME

def get_directories(dirs, exp_id, exp_no, tube_no):
    base_dir     = INI_DIR + '/'
    exp_base_dir = base_dir + exp_id + '/'
    exp_dir      = exp_base_dir + str(exp_no) + '/'
    tube_dir     = exp_dir + 'tube' + str(tube_no) + '/'
    masks_dir    = tube_dir + MASKS_DIRNAME + '/'

    if not dirs: dirs = ('search')

    single_result = False
    if type(dirs) == str:
         dirs = tuple(dirs)
         single_result = True
    results = []

    for dirtype in dirs:
        if dirtype == 'search':
            results.append((base_dir, exp_base_dir, exp_dir, tube_dir))
        elif dirtype == 'masks':
            results.append(masks_dir)
        elif dirtype == 'data':
            results.append(tube_dir)
        elif dirtype == 'base':
            results.append(base_dir)
        else:
            raise ValueError('Unknown value for get_directiories(): '
                             '{}'.format(dirs))

    if single_result: return results[0]
    else: return results

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
