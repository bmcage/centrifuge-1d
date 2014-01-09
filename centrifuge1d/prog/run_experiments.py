#!/usr/bin/python
from __future__ import print_function

import argparse
from os import listdir, sep
from ..shared import get_directories
from ..config import ModulesManager, load_model
from ..const import CSV_DIR, INI_DIR, FIGS_DIR

def compare_experiment_csv(exp_id, dataout_path=None):
    import os, subprocess
    from ..const import DIFFPROG

    if not DIFFPROG:
        print('No DIFFPROG was selected. Cannot compare. Exiting.')
        exit(1)

    exp_info = {key: '' for key in ['exp_no', 'mask']}
    exp_info['exp_id'] = exp_id
    if dataout_path:
        exp_info['ini_dir'] = dataout_path + sep + 'datafiles'
        exp_info['figs_dir'] = dataout_path + sep + 'datafiles' 
    else:
        exp_info['ini_dir'] = INI_DIR
        exp_info['figs_dir'] = FIGS_DIR

    csv_filename = ''
    exp_basedir = get_directories('ini', 'exp_base', exp_info)
    filenames = listdir(exp_basedir)
    for filename in filenames:
        if ((len(filename) > 10) and (filename[:5] == 'orig.')
            and (filename[-4:] == '.csv')):

            csv_filename = filename[5:]
            break

    if not csv_filename:
        print('Original csv file not found. Cannot compare. Exiting.')
        exit(0)

    csv_original = exp_basedir + '/' + 'orig.' + csv_filename
    csv_current = CSV_DIR + '/' + csv_filename

    if not os.path.exists(csv_current):
        csv_current = csv_current[:-4] + '/' + csv_filename

    if not os.path.exists(csv_current):
        print('Current CSV file: ', csv_current, 'does not exist. Cannot '
              'compare. Exiting.')
        exit(0)

    subprocess.call([DIFFPROG, csv_original, csv_current])

    exit(0)

def list_experiments(args, verbose=True, dataout_path=None):
    if not args.exp_id:
        list_type = 'base'
    elif not args.exp_no:
        list_type = 'exp_base'
    elif not args.mask:
    #     list_type = 'exp_no'
    # else:
        list_type = 'masks'

    exp_info = {key: getattr(args, key) for key in ['exp_id', 'exp_no', 'mask']}
    if dataout_path:
        exp_info['ini_dir'] = dataout_path + sep + 'datafiles'
        exp_info['figs_dir'] = dataout_path + sep + 'datafiles' 
    else:
        exp_info['ini_dir'] = INI_DIR
        exp_info['figs_dir'] = FIGS_DIR
    filenames = listdir(get_directories('ini', list_type, exp_info))

    if list_type == 'base':
        print('Available exp_IDs:\n')
        out_filenames = \
          [fname for fname in filenames if not fname[-4:] in ('.ini', '.csv')]

    elif list_type == 'exp_base':
        print("Available exp_NOs for experiment '" + args.exp_id)
        out_filenames = \
          [fname for fname in filenames if not fname[-4:] in ('.ini', '.csv')]

    # elif list_type == 'exp_no':
    #     print('Inifiles listed for experiment: ', args.exp_id,
    #           '\nwith exp_NO: ', args.exp_no)
    #     out_filenames = filter(lambda fname: not fname[:-4] == '.ini',
    #                            filenames)
    else:
        print("Available masks for experiment: '" + args.exp_id +
              "' of number '" + repr(args.exp_no) +"':\n")
        out_filenames = \
          [fname[:-4] for fname in filenames if fname[-4:] == '.ini']

    print('\n'.join(sorted(out_filenames)))

def list_modules(verbose=True):
    modman = ModulesManager()
    modman.echo(verbose)

def parse_input(dataout_path=None, parse_arg=[]):
    """
    Parse the calling sequence of the script. If run independent, 
    don't pass parse_arg and sys.argv will be used.
    Otherwise, pass in the arguments as a list, eg
     parse_arg = ['gem-zwijnaarde2-drain', '1', '-m', 'f_mo_limit']
    to run from folder gem-zwijnaarde2-drain, exp 1, with mast f_mo_limit.
    """
    usage_str = '\n\t%(prog)s [options] [exp_ID] [exp_NO] [last_exp_NO]'

    argparser = argparse.ArgumentParser(usage=usage_str)
    argparser.add_argument('exp_id', type=str, nargs="?",
                           metavar='exp_ID',
                           help=('ID identifying an experiment (see -l for all '
                                 'available experiments'))
    argparser.add_argument('exp_no', nargs="?", metavar='exp_NO',
                           help=("the number of the first experiment in the "
                                 "experiments serie identified by 'exp_ID' "))
    argparser.add_argument('last_exp_no', nargs='?', metavar='last_exp_NO',
                           help=("the number of the last experiment. All "
                                 "experiments between the 'exp_NO' "
                                 "and 'last_exp_NO' (included) will be "
                                 "performed. If not specified, only the "
                                 "experiment with 'exp_NO' is computed"))
    argparser.add_argument('-c', '--compare', dest='compare_p', default=False,
                           action='store_true',
                           help=("Compare two configuration files."))
    argparser.add_argument('-d', '--modules-list', dest='modules_list',
                           action="store_true", default=False,
                           help=("Get the list of all available centrifuge "
                                 "modules. Using the -v option gives also "
                                 "more detailed description of modules."))
    argparser.add_argument('--dry-run', dest='dry_run',
                           action="store_true", default=False,
                           help=("Simulate running, but without the actual "
                                 "computation. For debugging purposes only."))
    argparser.add_argument('-l', '--list', dest='list', action="store_true",
                           default=False,
                           help="List. If 'exp_ID' is not provided, list all "
                           "available experiments. If 'exp_NO' is not "
                           "provided, list all exp_NO under given 'exp_ID'. "
                           "If both provided, list all available mask files.")
    argparser.add_argument('-m', '--mask', dest='mask', action="append",
                           metavar='MASK_NAME', type=str,
                           help=("Run given experiment with a mask MASK_NAME "
                                 "file. The file \'MASK_NAME.ini\' in the masks "
                                 "directory should exist. Option can be supplied "
                                 "multiple times. In that case masks will be "
                                 "applied in the order they were specified"))
    argparser.add_argument('-p', '--print-config', dest='print_config_p',
                           action='store_true', default=False,
                           help=('Print the used configuration file for given '
                                 'experiment and exit.'))
    argparser.add_argument('-s', '--show', dest='show_p', default=False,
                           action='store_true',
                           help=("Show results if present. No computation is "
                                 "peformed."))
    argparser.add_argument('-v', '--verbose', dest='verbose',
                           action="store_true", default=False,
                           help=("Provide more detailed informations about "
                                 "what is being done."))
    argparser.add_argument('-w', '--csv-diff', dest='csv_diff', default=False,
                           action='store_true',
                           help=("Compare the original with actual csv file."))

    if parse_arg:
        args = argparser.parse_args(parse_arg)
    else:
        args = argparser.parse_args()

    if args.list:
        list_experiments(args, args.verbose)
        exit(0)
    elif args.modules_list:
        list_modules(args.verbose)
        exit(0)
    elif not args.exp_id:
        argparser.print_help()
        exit(0)
    elif args.csv_diff:
        compare_experiment_csv(args.exp_id, dataout_path)
        exit(0)
    elif args.exp_id and args.exp_no:
        pass # all needed data is passed
    else:
        argparser.print_help()
        exit(0)

    if not args.last_exp_no is None:
        try:
            int(args.exp_no)
            int(args.last_exp_no)
        except:
            print("\nERROR: When running a range of experiments, "
                  "first and last experiment numbers must be of type 'int'."
                  "\nCannot proceed, exiting...")
            exit(1)


    return args

def get_exp_no(options):
    if options.last_exp_no is None:
        yield options.exp_no
    else:
        first_experiment = int(options.exp_no)
        last_experiment  = int(options.last_exp_no)

        for exp_no in range(first_experiment, last_experiment+1):
            yield str(exp_no)

def run_experiments(options, dataout_path=None):
    (exp_id, mask) = (options.exp_id, options.mask)
    print_cfg_only = options.print_config_p
    verbose = not print_cfg_only

    if print_cfg_only:
        modman = None
    else:
        modman = ModulesManager()

    if verbose and options.last_exp_no:
        print('\n\n GENERAL experiments informations\n' + (60 * '-') +
              '\n Experiment ID   : ' + exp_id +
              '\n First experiment: ' + options.exp_no +
              '\n Last  experiment: ' + options.last_exp_no + '\n'+ (60 * '-'))

    for exp_no in get_exp_no(options):
        if verbose:
            header = ("Executing experiment {} number {}.".format(exp_id,
                                                                  exp_no))
            print('\n', len(header) * '=', '\n', header,
                  '\n', len(header) * '=')

        experiment_info =  {'exp_id': exp_id, 'exp_no': exp_no,
                            'mask': mask}
        if dataout_path:
            experiment_info['ini_dir'] = dataout_path + sep + 'datafiles'
            experiment_info['figs_dir'] = dataout_path + sep + 'datafiles' 
        else:
            experiment_info['ini_dir'] = INI_DIR
            experiment_info['figs_dir'] = FIGS_DIR

        model = load_model(experiment_info, display_only=print_cfg_only,
                           modman = modman)

        if print_cfg_only: continue

        solver_module = modman.find_module(model.exp_type, submodule='run')

        if options.dry_run:
            if not hasattr(solver_module, 'dry_run'):
                print("\nERROR: Module '" + model.exp_type + "' has not "
                      "specified 'dry_run' function. Exiting...\n")
                exit(1)

            solver_module.dry_run(model)
        else:
            solver_module.run(model)

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

def compare2configs(options, dataout_path=None):
    from ..config import load_configuration

    print('Add information for the second configuration file:')
    exp_id2  = input('Experiment ID  : ').strip()
    exp_no2  = input('Experiment No. : ').strip()
    if len(options.mask) > 1:
        mask2 = options.mask[1]
    else:
        mask2    = input('Mask (optional): ').strip()

    exp_info1 = {'exp_id': options.exp_id, 'exp_no': options.first_experiment,
                 'mask': options.mask[0]}
    exp_info2 = {'exp_id': exp_id2, 'exp_no': int(exp_no2), 'mask': mask2}
    if dataout_path:
        exp_info1['ini_dir'] = dataout_path + sep + 'datafiles'
        exp_info1['figs_dir'] = dataout_path + sep + 'datafiles'
        exp_info2['ini_dir'] = dataout_path + sep + 'datafiles'
        exp_info2['figs_dir'] = dataout_path + sep + 'datafiles'
    else:
        exp_info1['ini_dir'] = INI_DIR
        exp_info1['figs_dir'] = FIGS_DIR
        exp_info2['ini_dir'] = INI_DIR
        exp_info2['figs_dir'] = FIGS_DIR
        
    cfg1 = load_configuration(exp_info1)
    cfg2 = load_configuration(exp_info2)

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

def main(dataout_path=None, parse_arg=[]):
    options = parse_input(dataout_path, parse_arg)

    if options.show_p:
        from ..modules.shared.show import show_results

        experiment_info =  \
          {'exp_id': options.exp_id, 'exp_no': options.exp_no,
           'mask': options.mask}
        
        if dataout_path:
            experiment_info['ini_dir'] = dataout_path + sep + 'datafiles'
            experiment_info['figs_dir'] = dataout_path + sep + 'datafiles' 
        else:
            experiment_info['ini_dir'] = INI_DIR
            experiment_info['figs_dir'] = FIGS_DIR

        show_results(experiment_info)

    elif options.compare_p:
        compare2configs(options, dataout_path)
    else:
        run_experiments(options, dataout_path)    

if __name__ == "__main__":
    main()