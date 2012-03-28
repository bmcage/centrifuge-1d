#!/usr/bin/python
from sys import path as syspath, argv as sysargv
from os.path import exists
from common import load_modules, make_collector, print_by_tube
from config import read_cfgs, merge_flattened_cfgs, print_flattened_cfg
from config import ModelParameters
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
    optparser.add_option('-p', '--print-config', dest='print_config_p',
                         action='store_true', default=False,
                         help='Print the used configuration file for given'\
                         ' experiment and exit; if also parameter ''-t'' is'\
                         ' included, the config file for the tube is included'\
                         'too.')
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
        print('Error: Flag ''-l'' cannot be used with an argument.')
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

    default_ini = exp_inifiles_dir + '/defaults.ini'
    if not exists(default_ini):
        print('For experiments serie  %s the file ''defaults.ini'' does'
              ' not exist. '
              'Module name cannot be determined' % exp_id)
        exit(1)

    [default_cfg] = read_cfgs(default_ini, preserve_sections_p=False)

    #collector = make_collector(options.tubes)

    base_module = find_module('base')

    module = find_module(default_cfg['module'])


    if hasattr(module, 'generate_tubes_suffixes'):
        generate_tubes_suffixes = module.generate_tubes_suffixes
    else:
        generate_tubes_suffixes = base_module.generate_tubes_suffixes

    if not options.print_config_p:
        print('\n\n GENERAL experiments informations'
              '\n---------------------------------------------------------'
              '\n ID: %s'
              '\n First experiment: %s'
              '\n Last  experiment: %s '
              '\n---------------------------------------------------------'
              % (exp_id, first_experiment, last_experiment))

    (tubes_suffixes, tubes_identifiers) = generate_tubes_suffixes(options.tubes)

    for exp_no in range(first_experiment, last_experiment+1):
        group_default_ini = (exp_inifiles_dir + '/experiment_' + str(exp_no)
                             + '-defaults.ini')

        if exists(group_default_ini):
            [group_default_cfg] = read_cfgs(group_default_ini,
                                            preserve_sections_p = False)
        else:
            group_default_cfg = {}

        for (tube_suffix, tube_identifier) in zip(tubes_suffixes,
                                                  tubes_identifiers):
            if not options.print_config_p:
                print('\n===================================================='
                      '=============='
                      '\nExecuting experiment %s%s of the problem %s.'
                      % (str(exp_no), tube_identifier, exp_id))

            tube_default_ini = (exp_inifiles_dir + '/experiment_'
                                    + str(exp_no) + tube_suffix
                                    + '-defaults.ini')

            inifile = 'experiment_' + str(exp_no) + tube_suffix + '.ini'
            filename_ini = (exp_inifiles_dir + '/' + inifile)
            if not exists(filename_ini):
                print('File "%s" does not exist, skipping...' % inifile)
                continue

            if tube_suffix and exists(tube_default_ini):
                [tube_default_cfg] = read_cfgs(tube_default_ini,
                                               preserve_sections_p = False)
            else:
                tube_default_cfg = {}

            [filename_cfg] = read_cfgs(filename_ini, preserve_sections_p=False)


            cfg = merge_flattened_cfgs({}, default_cfg, group_default_cfg,
                                       tube_default_cfg, filename_cfg)

            if options.print_config_p:
                print_flattened_cfg(cfg)
                continue

            module.adjust_cfg(cfg)

            model = ModelParameters(cfg)

            results = module.solve(model)

            print('Results:\n', results)

        #collector('collect', data=results, tube_no=tube_no)

        #collector('print-by-tube', data=print_by_tube)

if __name__ == "__main__":
    (options, args) = parse_input()
    run_experiments(options, args)
