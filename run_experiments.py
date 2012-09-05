#!/usr/bin/python
from sys import path as syspath, argv as sysargv
from os import listdir
from os.path import exists
from shared import (make_collector, print_by_tube, get_directories, yn_prompt,
                    get_default_ini_filename)
from config import ModulesManager, ModelParameters, Configuration
from optparse import OptionParser
from const import DEFAULTS_ININAME, CONSTANTS_ININAME

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
                print('\n'.join(sorted(listdir(
                    get_directories('base', None, None, None, None)))))
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

def process_global_constants(cfg, consts_cfg):
    if not consts_cfg: return

    for (name, value) in consts_cfg.iterate_values():
        if not cfg.get_value(name):
            cfg.set_value(name, value)

def load_configuration(exp_id, exp_no, tube_no, mask=None):
    (search_dirs, data_dir, masks_dir) = \
      get_directories(['search', 'data', 'masks'], exp_id, exp_no, tube_no)

    filter_existing = lambda fnames: list(filter(lambda fname: exists(fname),
                                                 fnames))
    prefix_with_paths = lambda fname, dirs: map(lambda cfgdir: cfgdir + fname,
                                                dirs)

    defaults_files = filter_existing(prefix_with_paths(DEFAULTS_ININAME,
                                                       search_dirs))

    measurements_filenames = listdir(data_dir)
    measurements_files = []
    for fname in measurements_filenames:
        print(fname)
        # valid measurement files are *.ini (i.e. >4 chars filename)
        # except for 'defaults.ini'
        if ((fname == DEFAULTS_ININAME) or (len(fname) <= 4)
            or (fname[-4:] != '.ini')):
            continue

        measurements_files.append(data_dir + fname)

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

    # Handle CONSTANTS inifiles
    constants_files = filter_existing(prefix_with_paths(CONSTANTS_ININAME,
                                                        search_dirs))
    consts_cfg = None
    if constants_files:
        consts_cfg = Configuration().read_from_files(*constants_files)

    return (cfg, consts_cfg)

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
        collector = make_collector(options.tubes)

        modman = ModulesManager()

    for exp_no in range(first_experiment, last_experiment+1):
        for tube_no in tubes:
            if verbose:
                header = ("Executing experiment {} number {:d}, tube {}."
                          "".format(exp_id, exp_no, tube_no))
                print('\n', len(header) * '=', '\n', header,
                      '\n', len(header) * '=')

            (cfg, consts_cfg) = load_configuration(exp_id, exp_no,
                                                   tube_no, mask)

            if print_cfg_only:
                header = ("Configuration file of experiment '{}' number {:d}, "
                          "tube {}".format(exp_id, exp_no, tube_no))
                print("\n", header, '\n', len(header) * '-')
                cfg.echo()
                continue

            cfg.set_defaults(modman)

            # Assign global values not present in (or based on) configuration
            process_global_constants(cfg, consts_cfg)

            if not cfg.is_valid(modman):
                print('\n\nConfiguration is NOT VALID.\n'
                      'The full configuration is:')
                cfg.echo()
                exit(1)

            cfg.adjust_cfg(modman)

            solver_module = modman.find_module(cfg.get_value('exp_type'),
                                               submodule='run')

            model = ModelParameters(cfg)

            results = solver_module.run(model)

            collector('collect', data=results)

    if not print_cfg_only:
        print('Results summary:\n')
        print_fn = lambda x: print(x[0])
        collector('print', print_format_fn=print_fn)
        #collector('print-by-tube', data=print_by_tube)

if __name__ == "__main__":
    options = parse_input()

    print_cfg_only = options.print_config_p
    run_experiments(options.exp_id, options.first_experiment,
                    options.last_experiment, options.tubes, options.mask,
                    verbose=(not print_cfg_only), print_cfg_only=print_cfg_only)
