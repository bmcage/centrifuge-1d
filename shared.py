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

def get_directiories(dirs, exp_id, exp_no, tube_no):
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
        if dirstype == 'search':
            results.append(base_dir, exp_base_dir, exp_dir, tube_dir, masks_dir)
        elif dirtype == 'masks':
            results.append(masks_dir)
        elif dirtyp == 'data':
             results.append(tube_dir)
        else:
            raise ValueError('Unknown value for get_directiories(): '
                             '{}'.format(dirs))

    if single_result: return results[0]
    else return results

def load_configuration(exp_id, exp_no, tube_no, mask=None):
    (search_dirs, data_dir, masks_dir) = \
      get_directiories(['search', 'data', 'masks'], exp_id, exp_no, tube_no)

    filter_existing = lambda fnames: list(filter(lambda fname: exists(fname),
                                                 fnames))
    prefix_with_paths = lambda fname, dirs: map(lambda cfgdir: cfgdir + fname,
                                                dirs)
    defaults_files = filter_existing(prefix_with_paths(DEFAULTS_ININAME,
                                                       search_dirs))

    measurements_filenames = listdir(data_dir)
    masurements_files = []
    for fname in measurements_filenames:
        # valid measurement files are *.ini (i.e. >4 chars filename)
        # except for 'defaults.ini'
        if ((fname == DEFAULTS_ININAME) or (len(fname) <= 4)
            or (fname[-4:] != '.ini')):
            continue

        measurements_filenames.append(fname)

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

    cfg = Configuration().read_from_files(*cfg_files)

    # Handle CONSTANTS.ini files
    constants_files = filter_existing(prefix_with_paths('CONSTANTS.ini',
                                                        search_dirs))
    consts_cfg = None
    if constants_files:
        consts_cfg = Configuration().read_from_files(*constants_files)

    return (cfg, consts_cfg)

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
