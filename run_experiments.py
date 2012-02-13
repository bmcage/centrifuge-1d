#!/usr/bin/python
from sys import path as syspath, argv as sysargv

syspath.append('/'.join(['.', 'odes', 'build', 'lib.linux-x86_64-3.2']))

INIFILES_BASE_DIR = 'sources/inifiles'
TUBES_NUMBERS     = [1, 2, 4, 5]

def usage():
    print('\nUsage:'
          '\n  run_experiments <exp_ID> <first_experiment> [last_experiment]'
          '\n           or'
          '\n run_experiments -i <inifile>'
          '\n'
          '\n    exp_ID: experiment_ID - a unique identifier of experiment series'
          '\n    first_experiment: number of first experiment in exp_ID series'
          '\n             of experiments'
          '\n    last_experiment: if specified, computes all experiments between'
          '\n             the first_experiment and last_experiment (included);'
          '\n             if not specified, computes only first_experiment')

    exit(0)


def run_experiments(exp_id, first_experiment, last_experiment):
    from common import load_modules_names
    from config import read_cfgs, merge_cfgs, flatten_cfg
    from base import ModelParameters
    from os.path import exists

    exp_inifiles_dir = INIFILES_BASE_DIR + '/' + exp_id

    exp_ini_defaults = exp_inifiles_dir + '/defaults.ini'
    if exists(exp_ini_defaults):
        [default_cfg] = read_cfgs(ini_defaults, preserve_sections_p=True)
    else:
        default_cfg = ''

    find_module = load_modules_names('run')

    print('\n'
          '---------------------------------------------------------'
          '\n GENERAL experiments informations'
          '\nFirst experiment: %s'
          '\nLast  experiment: %s '
          % (first_experiment, last_experiment))

    for exp_no in range(first_experiment, last_experiment+1):
        print('\n==========================================================='
              '\nExecuting experiment %s of the %s problem.'
              % (exp_no, exp_id))
        for tube_no in TUBES_NUMBERS:
            print('\n   Executing experiment for tube number: ', tube_no)

            inifilename = (exp_inifiles_dir + '/experiment_' + str(exp_no)
                           + '-filter' + str(tube_no) +'.ini')

            [cfg] = read_cfgs(inifilename, preserve_sections_p=True)

            module = find_module(cfg['experiment-data']['exp_type'])

            model_cfg = merge_cfgs(module.base_cfg(), [default_cfg, cfg])
            model_cfg = flatten_cfg(model_cfg)

            module.adjust_cfg(model_cfg)
            if not module.check_cfg(model_cfg):
                raise ValueError('Check_cfg failed.')

            model = ModelParameters(model_cfg)

            module.solve(model)

        print('\nExperiment %s finished.' % exp_no)

def parse_input():

    arg_len = len(sysargv)

    if arg_len == 1 or (arg_len == 2 and (sysargv[1] in ['-h', '--help'])):
        usage()
    elif arg_len == 3 and sysargv[1] in ['-i', '--inifile']:
        # run_directly_solver(inifile = sysargv[2])
        raise NotImplemented
    elif arg_len == 3 or arg_len == 4:

        exp_id = sysargv[1]
        print(sysargv)

        try:
            first_experiment = int(sysargv[2])
            if arg_len == 4:
                last_experiment = int(sysargv[3])
            else:
                last_experiment = first_experiment
        except:
            raise ValueError('Input error: first and last experiment have to be'
                             ' integers. Wrong input: %s' % sysargv[1:])

        run_experiments(exp_id, first_experiment, last_experiment)
    else:
        print('\nWrong arguments for run_experiments: %s' % sysargv[1:])
        usage()
        exit(0)

if __name__ == "__main__":
    parse_input()


