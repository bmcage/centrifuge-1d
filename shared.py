import numpy as np
from const import INI_DIR, MASKS_DIRNAME, DEFAULTS_ININAME
from os import listdir
from os.path import exists
from config import Configuration

def yn_prompt(question_str):
    while True:
        answ = input(question_str)
        answ = answ.lower()
        if answ in ['', 'y', 'yes']: return True
        if answ in ['n', 'no']: return False
def get_experiment_base_dirs(exp_id):
    base_dir      = INI_DIR + '/'
    exp_base_dir  = base_dir + exp_id + '/'

    return (base_dir, exp_base_dir)

def get_experiment_dir(exp_base_dir, exp_no):
    return exp_base_dir + str(exp_no) + '/'

def get_tube_dirs(experiment_dir, tube_number):
    tube_dir  = experiment_dir + 'tube' + str(tube_number) + '/'
    masks_dir = tube_dir + MASKS_DIRNAME + '/'
    return (tube_dir, masks_dir)

def load_configuration_file(basedir, *inifilenames):
    existing_files = list(filter(lambda fname: exists(basedir + fname),
                                 inifilenames))
    if existing_files:
        return Configuration().read_from_files(basedir, *existing_files)
    else:
        return None

def load_cfg_defaults(*cfg_dirs):
    fullnames = list(map(lambda cfgdir: cfgdir + DEFAULTS_ININAME, cfg_dirs))
    return load_configuration_file('', *fullnames)

def load_measurements_cfg(meas_dir):
    filenames = listdir(meas_dir)

    for fname in (DEFAULTS_ININAME, MASKS_DIRNAME):
        if fname in filenames:
            filenames.remove(fname)

    return load_configuration_file(meas_dir, *filenames)

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
