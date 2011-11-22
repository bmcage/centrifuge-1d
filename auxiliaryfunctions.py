import const, getopt, sys, os
import numpy as np

def parse_centrifuge_input(args):
    """
    Read the ini file, use default if none given on command line
    """
    
    INIFILES = [const.INIFILE_DEFAULT, \
            const.INIFILE_DEFAULT_EXT]
    OUTPUTDIR = const.DATA_DIR
    SAVECONFIG = ''

    try:
        options, leftargs = getopt.getopt(args[1:],
                                          const.SHORTOPTS, const.LONGOPTS)
    except (getopt.GetoptError, msg):
        print(msg)
        # return without filling anything if we could not parse the args
        print("Error parsing the arguments: %s " % self.args[1:])
        sys.exit(0)
        if leftargs:
            print('run_centrifuge.py does not understand argument %s' % leftargs)
            sys.exit(0)
    
    for opt_ix in range(len(options)):
        option, value = options[opt_ix]
        if option in ( '-i', '--inifiles'):
            ini_input = value.split(',')
            if len(ini_input) == 1:
                INIFILES[0] = ini_input[0]
            else:
                INIFILES = ini_input
        elif option in ('-o', '--outputdir'):
            OUTPUTDIR = value
        elif option in ('-s', '--saveconfig'):
            SAVECONFIG = value

    return [INIFILES, OUTPUTDIR, SAVECONFIG]

def ensure_dir(f):
    d = os.path.dirname(f)
    if not os.path.exists(d):
        os.makedirs(d)

def save_data(filename, data_dict):
    ensure_dir(filename)
    np.savez(filename, **data_dict)

def load_data(filename):
    npzfile = np.load(filename)
    return npzfile
