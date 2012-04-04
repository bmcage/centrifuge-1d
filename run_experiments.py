#!/usr/bin/python
from sys import path as syspath, argv as sysargv
from os.path import exists
from shared import (make_collector, print_by_tube, generate_tubes_suffixes)
from config import ModulesManager, Configuration, ModelParameters
from optparse import OptionParser


syspath.append('/'.join(['.', 'odes', 'build', 'lib.linux-x86_64-3.2']))

INIFILES_BASE_DIR = 'sources/inifiles'
DEFAULT_TUBES     = [1, 2, 4, 5]

def parse_input():

    usage_str = \
      ('\n%prog [options] ID, first_experiment [last_experiment]'
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
                         help=("Run experiment only on selected tubes, default "
                               "is:\n %default"))
    optparser.add_option('-l', '--list', dest='list', action="store_true",
                         default=False,
                         help="Lists all available experiments")
    optparser.add_option('-p', '--print-config', dest='print_config_p',
                         action='store_true', default=False,
                         help=('Print the used configuration file for given '
                               'experiment and exit; if also parameter ''-t'' '
                               'is included, the config file for the tube is '
                               'included too.'))
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

def load_configuration(inifilename):
     if exists(inifilename):
        return Configuration().read_from_files(inifilename)
     else:
        return None

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
    default_cfg = load_configuration(default_ini)

    #collector = make_collector(options.tubes)

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
        group_default_cfg = load_configuration(group_default_ini)

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

            if tube_suffix:
                tube_default_cfg = load_configuration(tube_default_ini)
            else:
                tube_default_cfg = None

            experiment_cfg = load_configuration(filename_ini)

            cfg = Configuration().merge(default_cfg, group_default_cfg,
                                       tube_default_cfg, experiment_cfg)

            if options.print_config_p:
                cfg.echo()

            modman = ModulesManager()

            cfg.set_defaults(modman)

            if not cfg.is_valid(modman):
                print('\n\nConfiguration is NOT VALID.\n'
                      'The full configuration is:')
                cfg.echo()
                exit(1)

            cfg.adjust_cfg(modman)

            solver_module = modman.find_module(cfg.get_value('exp_type'),
                                               submodule='run')

            model = ModelParameters(cfg, modman)

            results = solver_module.solve(model)

            print('Results:\n', results)

        #collector('collect', data=results, tube_no=tube_no)

        #collector('print-by-tube', data=print_by_tube)

if __name__ == "__main__":
    (options, args) = parse_input()
    run_experiments(options, args)
