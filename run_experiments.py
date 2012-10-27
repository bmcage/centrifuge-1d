#!/usr/bin/python
from __future__ import print_function
from sys import path as syspath, argv as sysargv
from os import listdir
from os.path import exists
from shared import get_directories, yn_prompt
from config import ModulesManager, ModelParameters, Configuration
from optparse import OptionParser
from const import DEFAULTS_ININAME, CONSTANTS_ININAME

syspath.append('/'.join(['.', 'odes', 'build', 'lib.linux-x86_64-3.2']))

def list_experiments(verbose=True):
    from os import listdir

    print('\n'.join(sorted(listdir(
        get_directories('ini', 'base',
                        {key: '' for key in ['exp_id', 'exp_no',
                                             'mask']})))))

def list_modules(verbose=True):
    if modules_list:
        modman = ModulesManager()
        modman.echo(verbose)

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
    optparser.add_option('-c', '--compare', dest='compare_p', default=False,
                         action='store_true',
                         help=("Compare two configuration files."))
    optparser.add_option('-d', '--modules-list', dest='modules_list',
                         action="store_true", default=False,
                         help=("Get the list of all available centrifuge "
                               "modules"))
    optparser.add_option('-p', '--print-config', dest='print_config_p',
                         action='store_true', default=False,
                         help=('Print the used configuration file for given '
                               'experiment and exit.'))
    optparser.add_option('-s', '--show', dest='show_p', default=False,
                         action='store_true',
                         help=("Show results if present. No computation is "
                               "peformed."))
    optparser.add_option('-v', '--verbose', dest='verbose',
                         action="store_true", default=False,
                         help="If possible, provide more detailed informations")

    (options, args) = optparser.parse_args()
    arg_len = len(args)
    if arg_len == 0:
        if options.list: or :
            list_experiments(options.verbose)
        elif options.modules_list:
            list_modules(options.verbose)
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

    except:

        raise ValueError('Input error: first and last experiment have to be'
                         ' integers. Wrong input: %s '
                         % args[1:])

    return options

def process_global_constants(cfg, consts_cfg):
    if not consts_cfg: return

    tube_no = cfg.get_value('tube_no')
    cfg.set_parameters(consts_cfg.get_value('tubes')[tube_no])

def load_configuration(experiment_info):
    (search_dirs, data_dir, masks_dir) = \
      get_directories('ini', ['search', 'data', 'masks'], experiment_info)

    filter_existing = lambda fnames: list(filter(lambda fname: exists(fname),
                                                 fnames))
    prefix_with_paths = lambda fname, dirs: map(lambda cfgdir: cfgdir + fname,
                                                dirs)

    defaults_files = filter_existing(prefix_with_paths(DEFAULTS_ININAME,
                                                       search_dirs))

    measurements_filenames = listdir(data_dir)
    measurements_files = []
    for fname in measurements_filenames:
        # valid measurement files are *.ini (i.e. >4 chars filename)
        # except for 'defaults.ini'
        if ((fname == DEFAULTS_ININAME) or (len(fname) <= 4)
            or (fname[-4:] != '.ini')):
            continue

        measurements_files.append(data_dir + fname)

    mask_filename = ''
    mask = experiment_info['mask']
    if mask:
        mask_filename = masks_dir + mask + '.ini'
        if not exists(mask_filename):
            print('Mask file "{}" does not exist in expected location:'
                  '\n{}.'.format(mask, masks_dir))
            if not yn_prompt('Do you wish to continue without applying '
                         'the mask? [Y/n]: '):
                exit(0)

            mask_filename = ''

    cfg_files = defaults_files + measurements_files + [mask_filename]

    cfg = Configuration().read_from_files(*cfg_files)

    # Handle CONSTANTS inifiles
    constants_files = filter_existing(prefix_with_paths(CONSTANTS_ININAME,
                                                        search_dirs))
    consts_cfg = None
    if constants_files:
        consts_cfg = Configuration().read_from_files(*constants_files)

    return (cfg, consts_cfg)

def run_experiments(exp_id, first_experiment, last_experiment, mask,
                    verbose=True, print_cfg_only=False):

    if verbose:
        print('\n\n GENERAL experiments informations'
              '\n---------------------------------------------------------'
              '\n ID              : %s'
              '\n First experiment: %s'
              '\n Last  experiment: %s'
              '\n---------------------------------------------------------'
              % (exp_id, first_experiment, last_experiment))

    if not print_cfg_only:
        modman = ModulesManager()

    for exp_no in range(first_experiment, last_experiment+1):
            if verbose:
                header = ("Executing experiment {} number {:d}."
                          "".format(exp_id, exp_no))
                print('\n', len(header) * '=', '\n', header,
                      '\n', len(header) * '=')

            experiment_info =  {'exp_id': exp_id, 'exp_no': exp_no,
                                'mask': mask}

            (cfg, consts_cfg) = load_configuration(experiment_info)

            if print_cfg_only:
                header = ("Configuration file of experiment '{}' number {:d}"
                          .format(exp_id, exp_no))
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

            model.experiment_info = experiment_info

            results = solver_module.run(model)

def iterate_value(arg):

    if type(arg) in [list, tuple]:
        for value in arg:
            yield value
    elif type(arg) == dict:
        for (key, value) in arg.items():
            # dirty hack: we assume that only dict of list value is
            # 'inv_init_params', so the value is (init, [range])
            if type(value) in [tuple, list]:
                yield key + ': ' + str(value[0])
            else:
                yield key + ': ' + str(value)
    else:
        yield arg

    while True:
        yield ''

def compare2configs(options):
    print('Add information for the second configuration file:')
    exp_id2  = input('Experiment ID  : ').strip()
    exp_no2  = input('Experiment No. : ').strip()
    mask2    = input('Mask (optional): ').strip()

    exp_info1 ={'exp_id': options.exp_id, 'exp_no': options.first_experiment,
                'mask': options.mask}
    exp_info2 ={'exp_id': exp_id2, 'exp_no': int(exp_no2), 'mask': mask2}

    (cfg1, const_cfg1) = load_configuration(experiment_info1)
    (cfg2, const_cfg2) = load_configuration(experiment_info2)

    all_parameters = sorted(set(cfg1.list_options() + cfg2.list_options()))

    param_max_length = 15
    value_max_length = 20
    fstring = (' {:' + str(param_max_length) +'} | {:'
               + str(value_max_length) +'} | {:'
               + str(value_max_length) + '}')
    print('\n', (param_max_length + 2*value_max_length + 6) * '-')
    print(fstring.format('Option name', '1. config', '2. config'),
          '\n', (param_max_length + 2*value_max_length + 6) * '-')

    for parameter in all_parameters:
        v1 = cfg1.get_value(parameter, not_found='')
        v2 = cfg2.get_value(parameter, not_found='')

        if not v1 == v2:
            v1_iter = iterate_value(v1)
            v2_iter = iterate_value(v2)

            if len(parameter) > param_max_length:
                parameter = parameter[:param_max_length]
            print(fstring.format(parameter, next(v1_iter), next(v2_iter)))

            while True:
                v1_value = next(v1_iter)
                v2_value = next(v2_iter)
                if (v1_value == '') and (v2_value == ''): break

                print(fstring.format('', v1_value, v2_value))



if __name__ == "__main__":
    options = parse_input()

    if options.show_p:
        from modules.shared.show import ResultsData, DPlots

        experiment_info =  \
          {'exp_id': options.exp_id, 'exp_no': options.first_experiment,
           'mask': options.mask}
        data = ResultsData()
        if data.load(experiment_info):
            dplots = DPlots(data, experiment_info)
            dplots.display()

    elif options.compare_p:
        compare2configs(options)
    else:
        print_cfg_only = options.print_config_p
        run_experiments(options.exp_id, options.first_experiment,
                        options.last_experiment, options.mask,
                        verbose=(not print_cfg_only),
                        print_cfg_only=print_cfg_only)
