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

def load_configuration_file(basedir, *inifilenames):
    existing_files = list(filter(lambda fname: exists(basedir + fname),
                                 inifilenames))
    if existing_files:
        return Configuration().read_from_files(basedir, *existing_files)
    else:
        return None

def load_configuration(exp_id, exp_no, tube_no, mask=None):
    base_dir     = INI_DIR + '/'
    exp_base_dir = base_dir + exp_id + '/'
    exp_dir      = exp_base_dir + str(exp_no) + '/'
    tube_dir     = exp_dir + 'tube' + str(tube_no) + '/'
    masks_dir    = tube_dir + MASKS_DIRNAME + '/'

    search_dirs  = (base_dir, exp_base_dir, exp_dir, tube_dir)

    defaults_files = list(filter(lambda fname: exists(fname),
                                 map(lambda cfgdir: cfgdir + DEFAULTS_ININAME,
                                     search_dirs)))

    # Measurements: 1. list filenames 2. remove 'defaults.ini' and 'mask' dir
    #               3. prefix filenames with path
    measurements_filenames = listdir(tube_dir)
    for fname in (DEFAULTS_ININAME, MASKS_DIRNAME):
        if fname in measurements_filenames:
            measurements_filenames.remove(fname)
    measurements_files = list(map(lambda fname: tube_dir + fname,
                                  measurements_filenames))

    if mask:
        mask_filename = masks_dir + mask + '.ini'
        if not exists(mask_filename):
            print('Mask file "{}" does not exist in expected location:'
                  '\n{}.'.format(mask, masks_dir))
            if yn_prompt('Do you wish to continue without applying '
                         'the mask? [Y/n]: '):
                mask_filename = None
            else:
                exit(0)

    cfg_files = defaults_files + measurements_files + [mask_filename]

    return Configuration().read_from_files(*cfg_files)

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
