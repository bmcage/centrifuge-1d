from os.path import exists

def experiments_files(first_experiment, last_experiment, tubes):
    files = []
    identifiers = []

    for exp_no in range(first_experiment, last_experiment+1):
        for tube_no in tubes:
            inifilename = ('experiment_' + str(exp_no)
                           + '-filter' + str(tube_no) +'.ini')
            identifier  = 'experiment ' + str(exp_no) + ', tube ' + str(tube_no)

            files.append(inifilename)
            identifiers.append(identifier)

    return (identifiers, files)
