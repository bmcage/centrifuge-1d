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
    from common import modules_names

    print('\n'
          '---------------------------------------------------------'
          '\n GENERAL experiments informations'
          '\n'
          'First experiment: %s'
          'Last  experiment: %s '
          % first_experiment, last_experiment)

    inifiles_dir = ''.join(INIFILES_BASE_DIR, '/', exp_id)

    find_module = modules_names()

    for exp_no in range(first_experiment, last_experiment+1):
        for tube_no in TUBES_NUMBERS:
            filename = ''.join(inifiles_dir, '/experiment_', tube_no, '.ini')

            print('\n============================================================='
                  '\nExecuting experiment %s of the %s problem.'
                  % tube_no, exp_id)

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

        first_experiment = sysargv[3]
        if arg_len == 4:
            last_experiment = sysargv[4]
        else:
            last_experiment = first_experiment

        run_experiments(exp_id, first_experiment, last_experiment)
    else:
        print('\nWrong arguments for run_experiments: %s' % sysargv[1:])
        usage()
        exit(0)

if __name__ == "__main__":
    parse_input()


