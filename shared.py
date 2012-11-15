from __future__ import print_function
from const import INI_DIR, MASKS_DIRNAME, FIGS_DIR

def yn_prompt(question_str, default='y'):
    while True:
        try:
            answ = raw_input(question_str)
        except:
            answ = input(question_str)
        if answ == '': answ = default
        else: answ = answ.lower()
        if answ in ['y', 'yes']: return True
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

################################################################################
#                            Parsing functions                                 #
################################################################################

def parse_list(raw_value):
    if raw_value == '': return []

    result = []
    brackets = ('(', ')', '[', ']', '{', '}')
    counts = {br: 0 for br in brackets}

    i0 = 0
    for (i, val) in enumerate(raw_value):
        if val == ',':
            if ((counts['('] == counts[')']) and (counts['{'] == counts['}'])
                and (counts['['] == counts[']'])): # balanced brackets

                result.append(parse_value(raw_value[i0:i]))
                i0 = i+1
            else:
                pass
        elif val in brackets:
            counts[val] += 1

    result.append(parse_value(raw_value[i0:])) # append last value

    return result

def parse_dict(raw_value):
    if raw_value == '': return {}

    result = {}
    brackets = ('(', ')', '[', ']', '{', '}')
    counts = {br: 0 for br in brackets}

    i0 = 0
    for (i, char) in enumerate(raw_value):
        if char == ',':
            if ((counts['('] == counts[')']) and (counts['{'] == counts['}'])
                and (counts['['] == counts[']'])): # balanced brackets

                k = i0 + raw_value[i0:].index(':')
                if k>i:
                    print('Could not parse dict item, no key:value',
                          'value pair found:\n', raw_value[i0:i],
                          '\nof item:\n', raw_value)

                key   = parse_value(raw_value[i0:k])
                value = parse_value(raw_value[k+1:i])
                result[key] = value

                i0 = i+1
            else:
                pass
        elif char in brackets:
            counts[char] += 1

    k     = i0 + raw_value[i0:].index(':')
    key   = parse_value(raw_value[i0:k])
    value = parse_value(raw_value[k+1:])
    result[key] = value

    return result

def parse_value(str_value):
    """
      Given a value as string, tries to convert to it's correspending type.
      May be called recursively in the case of nested structures.
    """
    try:
        raw_value = str_value.strip()

        if raw_value[0] == "[" and raw_value[-1] == "]":
            return parse_list(raw_value[1:-1])
        elif raw_value[0] == "(" and raw_value[-1] == ")":
            return tuple(parse_list(raw_value[1:-1]))
        elif raw_value[0] == "{" and raw_value[-1] == "}":
            return parse_dict(raw_value[1:-1])
        elif ((raw_value[0] == "'" and raw_value[-1] == "'")
              or (raw_value[0] == '"' and raw_value[-1] == '"')):
            return raw_value[1:-1]
        elif raw_value == 'True':
            return True
        elif raw_value == 'False':
            return False
        elif raw_value == 'None':
            return None
        elif raw_value == 'inf':
            return inf
        elif raw_value == '-inf':
            return -inf
        elif '*' in raw_value:
            [raw_mul, raw_val] = raw_value.split("*")
            mul = parse_value(raw_mul)
            val = parse_value(raw_val)
            return [val for i in range(mul)]
        elif ':' in raw_value:
            range_values = raw_value.split(':')
            rstart = parse_value(range_values[0])
            rstop  = parse_value(range_values[1]) + 1
            if len(range_values) > 2:
                rstep = parse_value(range_values[2])
            else:
                rstep = 1
            return list(range(rstart, rstop, rstep))
        elif "." in raw_value or "e" in raw_value or "E" in raw_value:
            return float(raw_value)
        else:
            return int(raw_value)

    except:
        print('Error:Could not parse value: ', str_value, '\nExiting...')
        exit(1)
