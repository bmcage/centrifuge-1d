#!/usr/bin/python
from sys import path as syspath, argv as sysargv
from os.path import exists
from common import load_modules, make_collector, print_by_tube
from config import read_cfgs, merge_cfgs, flatten_cfg, ModelParameters
from optparse import OptionParser


syspath.append('/'.join(['.', 'odes', 'build', 'lib.linux-x86_64-3.2']))

INIFILES_BASE_DIR = 'sources/inifiles'
DEFAULT_TUBES     = [1, 2, 4, 5]

find_module = load_modules('run')

def parse_input():

    usage_str = ('\n%prog [options] ID, first_experiment [last_experiment]'
             '\n\n\tfirst_experiment:'
             '\n\t\tnumber of the first experiment in exp_ID series of experiments;'
             '\n\t\teither "first_experiment" or "-i" option has to be set'
             '\n\tlast_experiment:'
             '\n\t\tif specified, computes all experiments between the'
             '\n\t\t"first_experiment" and the "last_experiment" (included);'
             '\n\t\tif not specified, computes only first_experiment')
    optparser = OptionParser(usage=usage_str)
    optparser.add_option('-t', '--tubes', dest='tubes', default=DEFAULT_TUBES,
                         metavar='TUBES_NUMBERS',
                         help="Run experiment only on selected tubes, default is:\n %default")
    optparser.add_option('-l', '--list', dest='list', action="store_true",
                         default=False,
                         help="Lists all available experiments")
    (options, args) = optparser.parse_args()
    arg_len = len(args)
    if arg_len == 0:
        if options.list:
            from os import listdir
            print('\n'.join(listdir(INIFILES_BASE_DIR)))
            exit(0)
        optparser.print_help()
        exit(0)
    elif arg_len == 1:
        optparser.print_help()
        exit(0)
    elif options.list:
        print('Error: Flag ''-l'' cannot be used with and argument.')
        exit(1)

    optparser.destroy()

    return (options, args)

def run_experiments(options, exp_args):

    try:
        exp_id = exp_args[0]
        first_experiment = int(exp_args[1])
        if len(exp_args) == 3:
            last_experiment = int(exp_args[2])
        else:
            last_experiment = first_experiment
    except:

        raise ValueError('Input error: first and last experiment have to be'
            ' integers. Wrong input: %s ' 
            % exp_args[1:])

    exp_inifiles_dir = INIFILES_BASE_DIR + '/' + exp_id

    ini_defaults = exp_inifiles_dir + '/defaults.ini'
    if not exists(ini_defaults):
        print('Experiment %s does not have file defaults.ini'
              'Module name cannot be determined' % exp_id)
        exit(0)

    [default_cfg] = read_cfgs(ini_defaults, preserve_sections_p=True)

    #collector = make_collector(options.tubes)

    default_module = find_module('base')

    module = find_module(default_cfg['solver']['module'])

    if hasattr(module, 'experiments_files'):
        experiments_files = module.experiments_files
    else:
        experiments_files = default_module.experiments_files

    print('\n\n GENERAL experiments informations'
        '\n---------------------------------------------------------'
        '\n ID: %s'
        '\n First experiment: %s'
        '\n Last  experiment: %s '
        '\n---------------------------------------------------------'
        % (exp_id, first_experiment, last_experiment))

    (identifiers, experiments_inifiles) = \
      experiments_files(first_experiment, last_experiment, options.tubes)

    for (identifier, inifile) in zip(identifiers, experiments_inifiles):
        print('\n=========================================================='
            '\nExecuting %s of the problem %s.' % (identifier, exp_id))

        inifile_fullname = exp_inifiles_dir + '/' + inifile
        if not exists(inifile_fullname):
            print('File "%s" does not exist, skipping...' % inifile)
            continue

        [cfg] = read_cfgs(inifile_fullname, preserve_sections_p=True)

        model_cfg = merge_cfgs(module.base_cfg(), [default_cfg, cfg])
        model_cfg = flatten_cfg(model_cfg)

        module.adjust_cfg(model_cfg)

        model = ModelParameters(model_cfg)

        results = module.solve(model)

        print(results)

        #collector('collect', data=results, tube_no=tube_no)

        #collector('print-by-tube', data=print_by_tube)

if __name__ == "__main__":
    (options, args) = parse_input()
    run_experiments(options, args)


