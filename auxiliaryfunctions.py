import os
import numpy as np

def ensure_dir(f):
    d = os.path.dirname(f)
    if not os.path.exists(d):
        os.makedirs(d)

def ensure_numpy_array(value, arraysize_if_single_value = 1):
    try:
        if np.isscalar(value):
            array    = np.empty([arraysize_if_single_value, ], float)
            array[:] = value

            return array
        else:
            return np.asarray(value, float)
    except:
        raise ValueError('af.ensure_numpy_value: value is not a number or '
                         'a sequence of numbers: %s' % value)

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
