import numpy as np

def load_modules(submodule_name):
    """
    Loads all submodules of name 'submodule_name', which are accessed by
    the module name. Loading is done only on demand.
    """

    from os import listdir
    from sys import modules as sysmodules

    modules_names = listdir('modules')
    modules_names.remove('__pycache__')
    modules_names.remove('__init__.py')
    modules_names.remove('shared')

    loaded_modules = dict.fromkeys(modules_names)

    def find_module(module_name):
        if not module_name in loaded_modules:
            print('Module load error: Module "%s" could not be loaded.'
                      ' Unknown module name.\n'  % module_name)
            exit(1)

        module = loaded_modules[module_name]

        if not module:
            module_full_name = 'modules.' + module_name + '.' + submodule_name
            __import__(module_full_name)

            module = sysmodules[module_full_name]
            loaded_modules[module_name] = module

        return module

    return find_module

def load_experiment_types():

    from os import listdir
    from sys import modules as sysmodules

    modules_names = listdir('modules')
    modules_names.remove('__pycache__')
    modules_names.remove('__init__.py')
    modules_names.remove('shared')

    types_to_names = {}

    for module_name in modules_names:
        try:
            module_info_full_name = 'modules.' + module_name + '.info'

            __import__(module_info_full_name)

            info = sysmodules[module_info_full_name]

            for exp_type in info.types:
                types_to_names[exp_type] = module_name
        except:
            print('Module load error:Module "%s" could not be loaded.'
                  ' Skipping.\n'  % module_name)

    def find_module_name(exp_type):
        if not exp_type in types_to_names:
            print('Experiment error: Unknown experiment type: "%s"'
                  % module_name)
            exit(1)

        return types_to_names[exp_type]

    return find_module_name

def generate_tubes_suffixes(tubes_numbers):
    if not tubes_numbers:
        print('generate_tubes_suffixes error: No tubes (numbers) were '
              'specified - input values was: "%s"' % tubes_numbers)
        exit(1)

    suffixes = ['-tube' + str(tube_no) for tube_no in tubes_numbers]
    identifiers = [', tube ' + str(tube_no) for tube_no in tubes_numbers]

    return suffixes, identifiers


def print_by_tube(tube_number, tube_data):
    print('Tube number: ', tube_number)
    for tdata in tube_data:
        print(tdata[0] / 100.)

def make_collector(tubes_numbers):
    by_tube_collection = {tube_no:[] for tube_no in tubes_numbers}
    fifo_collection = []

    def collection(command, data = None, tube_no = None):
        if command == 'collect':
            if tube_no:
                by_tube_collection[tube_no].append(data)
            else:
                raise ValueError('Collector:collect - tube number or data'
                                 ' not specified')
        elif command in ['print', 'print-by-tube', 'print-fifo']:
            if data:
                if command == 'print-by-tube':
                    for (tube_no, tube_data) in by_tube_collection.items():
                        data(tube_no, tube_data)
                else:
                    data(fifo_collection)
            else:
                raise ValueError('Collector:print data not specified. Data'
                                 ' has to contain a print-data function.')
        elif command == 'get':
            return (fifo_collection, by_tube_collection)
        else:
            raise ValueError('Collector: Unknown command: %s.\n Valid commans'
                             ' are: "collect", "print", "print-by-tube",'
                             ' "print-fifo" or "get".' % command)
    return collection
