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

def ensure_numpy_array(value):
    try:
        if (type(value) == float or type(value) == int
            or type(value) == np.float
            or type(value) == np.float32
            or type(value) == np.float64
            or type(value) == np.float128
            or type(value) == np.int
            or type(value) == np.int8
            or type(value) == np.int16
            or type(value) == np.int32
            or type(value) == np.int64):

            return np.asarray([value], float)
        else:
            return np.asarray(value, float)
    except:
        raise ValueError('af.ensure_numpy_value: value not a number or sequence of numbers: %s' % value)

def save_data_measured_saturated_heights(filename, t, h1, h2, L, omega, porosity):
    t  = ensure_numpy_array(t)
    h1 = ensure_numpy_array(h1)
    h2 = ensure_numpy_array(h2)
    L  = ensure_numpy_array(L)
    omega = ensure_numpy_array(omega)
    porosity = ensure_numpy_array(porosity)

    save_data(filename, {'t': t,
                         'heights0': h1, 'heights1': h2,
                         'length': L,
                         'omega': omega,
                         'porosity': porosity})

def load_data(filename):
    npzfile = np.load(filename)
    return npzfile

def apply_functions(fns, *args):
    """
    Calls functions in 'fns', each with argument(s) given in 'args'.
    Functions 'fns' are assumed to be a single function or a LIST of functions
    and that the arguments 'args' can be modified in place(as no output of fns
    is stored and the function 'apply_functions' has no return value.

    Functions are evaluated in reversed order, i.e. the last function in 'fns'
    is applied first and the first in list is applied as the las one.

    See also: compose_functions
    """
    if fns:
        if type(fns) == list:
            fns.reverse()
            for fn in fns:
                fn(*args)
            fns.reverse()
        else:
            fns(*args)

def compose_functions(fns, data):
    """
    Calls functions in 'fns', each with argument(s) given in 'args'.
    Functions 'fns' are assumed to be a single function or a LIST of functions.
    The next function is called with argument, that is the result of the call
    of the previous function. The result is returned.

    Functions are evaluated in reversed order, i.e. the last function in 'fns'
    is applied first and the first in list is applied as the las one.

    See also: apply_functions
    """
    if fns:
        if type(fns) == list:
            result = data
            fns.reverse()
            for fn in fns:
                result = fn(result)
            fns.reverse()
        else:
            fns(data)

        return result

def load_centrifuge_configs(inifilenames, post_hook_fns = None):
    """
    Loads centrifuge inifiles (filenames - either a file or a list of files)
    and on the loaded data object calls the post_hook_fn (if set).
    """
    from config import ConfigManager, DEFAULT_PARAMETERS

    cm = ConfigManager.get_instance()
    model   = cm.read_configs(merge_output_p = True,
                              preserve_sections_p = False,
                              filenames = inifilenames,
                              defaults = [DEFAULT_PARAMETERS],
                              saveconfigname = '')
    apply_functions(post_hook_fns, model)

    return model

def load_measured_data(filename, post_hook_fns = None):
    from config import ConfigManager, DEFAULT_DATA_PARAMETERS

    cm = ConfigManager.get_instance()
    measured_data = cm.read_configs(merge_output_p = True,
                                    preserve_sections_p = False,
                                    filenames = [filename],
                                    defaults = [DEFAULT_DATA_PARAMETERS],
                                    saveconfigname = '')
    apply_functions(post_hook_fns, measured_data)

    return measured_data
