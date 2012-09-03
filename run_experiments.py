#!/usr/bin/python
from sys import path as syspath, argv as sysargv
from shared import (make_collector, print_by_tube, load_configuration)
from config import ModulesManager, ModelParameters
from optparse import OptionParser

syspath.append('/'.join(['.', 'odes', 'build', 'lib.linux-x86_64-3.2']))

DEFAULT_TUBES     = '1,2,4,5'

def parse_input():

    usage_str = \
      ('\n  %prog [options] exp_ID, first_experiment[, last_experiment]'
       '\n\nArguments:'
       '\n  exp_ID:'
       '\n      ID identifying an experiment (see -l for all available '
       'experiments)'
       '\n  first_experiment:'
       '\n      the number of the first experiment in the \'exp_ID\' '
       'experiments serie'
       '\n  last_experiment:'
       '\n      all experiments between the \'first_experiment\' and '
       'the \'last_experiment\' (included) will be computed.'
       '\n      if it is not specified, computes only the \'first_experiment\'')
    optparser = OptionParser(usage=usage_str)
    optparser.add_option('-l', '--list', dest='list', action="store_true",
                         default=False,
                         help="Lists all available experiments")
    optparser.add_option('-m', '--mask', dest='mask', default='',
                         metavar='MASK_NAME',
                         help=("Run given experiment with a mask MASK_NAME "
                               "file. Used in conjuction with ''-t'' parameter."))
    optparser.add_option('-d', '--modules-list', dest='modules_list',
                         action="store_true", default=False,
                         help=("Get the list of all available centrifuge "
                               "modules"))
    optparser.add_option('-p', '--print-config', dest='print_config_p',
                         action='store_true', default=False,
                         help=('Print the used configuration file for given '
                               'experiment and exit; if also parameter ''-t'' '
                               'is included, the config file for the tube is '
                               'included too.'))
    optparser.add_option('-t', '--tubes', dest='tubes', default=DEFAULT_TUBES,
                         metavar='TUBES_NUMBERS',
                         help=("Run experiment only on selected tubes, default "
                               "is:\n %default"))
    optparser.add_option('-v', '--verbose', dest='verbose',
                         action="store_true", default=False,
                         help="If possible, provide more detailed informations")

    (options, args) = optparser.parse_args()
    arg_len = len(args)
    if arg_len == 0:
        if options.list or options.modules_list:
            from os import listdir

            if options.list:
                print('\n'.join(sorted(listdir(INI_DIR))))
            if options.modules_list:
                modman = ModulesManager()
                modman.echo(options.verbose)
        else:
            optparser.print_help()

        exit(0)
    elif arg_len == 1:
        optparser.print_help()
        exit(0)
    elif options.list:
        print('Error: Flag ''-l'' cannot be used with an argument.')
        exit(1)

    optparser.destroy()

    try:
        exp_id = args[0]
        first_experiment = int(args[1])
        if len(args) == 3:
            last_experiment = int(args[2])
        else:
            last_experiment = first_experiment

        options.exp_id           = exp_id
        options.first_experiment = first_experiment
        options.last_experiment  = last_experiment
        options.tubes            = options.tubes.split(',')

    except:

        raise ValueError('Input error: first and last experiment have to be'
                         ' integers. Wrong input: %s '
                         % args[1:])

    return options

def run_experiments(exp_id, first_experiment, last_experiment, tubes, mask,
                    verbose=True, print_cfg_only=False):

    if verbose:
        print('\n\n GENERAL experiments informations'
              '\n---------------------------------------------------------'
              '\n ID              : %s'
              '\n First experiment: %s'
              '\n Last  experiment: %s'
              '\n Tubes           : %s'
              '\n---------------------------------------------------------'
              % (exp_id, first_experiment, last_experiment, ','.join(tubes)))

    if not print_cfg_only:

        modman = ModulesManager()

    for exp_no in range(first_experiment, last_experiment+1):

        for tube_no in tubes:
            if verbose:
                print('\n', 60 * '=',
                      '\n   Executing experiment %s%s of the problem ''%s''.'
                      % (str(exp_no), 'tube %s' % tube_number, exp_id),
                      '\n', 60 * '=')


            cfg = load_configuration(exp_id, exp_no, tube_no, mask)


            cfg.set_defaults(modman)


            cfg.adjust_cfg(modman)




            solver_module = modman.find_module(cfg.get_value('exp_type'),
                                               submodule='run')

            model = ModelParameters(cfg)

            results = solver_module.run(model)

            collector('collect', data=results)


if __name__ == "__main__":
    options = parse_input()

    print_cfg_only = options.print_config_p
    run_experiments(options.exp_id, options.first_experiment,
                    options.last_experiment, options.tubes, options.mask,
                    verbose=(not print_cfg_only), print_cfg_only=print_cfg_only)
