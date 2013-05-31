from __future__ import division, print_function

import numpy as np
import matplotlib.pyplot as plt
import functions

class Data():
    def __init__(self):
        pass

def main():
    #read from argument the config file to show
    import sys
    import os.path
    try:
        import ConfigParser as configparser
    except:
        import configparser
    ARGS = sys.argv
    if len(ARGS) < 2 or ARGS[1] in ['-h', '--help', 'help']:
        print ("""
USAGE: Give as first argument ini file with gramforce measurement data.
       After this, optionally give couples of options with name, like:
            * xmin value
            * xmax value

Example use:
$ python prog_plot_meas.py ~/git/centrifuge-1d/data/datafiles/gem-mixture-drain/16/Data\ Instr\ INSTR\ 9_27_2012\ 14_07_51.ini xmin 800 xmax 1000

""")
        sys.exit(0)
    fname = ARGS[1]
    if not os.path.isfile(fname):
        print ("ERROR: %s is not a file on your PC" % fname)
        sys.exit(0)
    optname = []
    optval = []
    for opt in ARGS[2::2]:
        optname.append(opt)
    for opt in ARGS[3::2]:
        optval.append(opt)

    #now we open the file
    parser   = configparser.ConfigParser()
    if sys.version_info[0] < 3:
        try:
            parser.read(fname)
        except Exception as E:
            print(E)
            exit(0)
    else:
        try:
            parser.read(fname)
        except configparser.DuplicateOptionError as E:
            print(E)
            exit(0)

    # Write data from parser to configuration
    data = Data()
    for psection in parser.sections():
        section = psection.lower()

        for option in parser.options(psection):
            raw_value = parser.get(psection, option).strip()
            setattr(data, option.lower(), raw_value)
    ind = 0
    for prefix, label in [('gf_mo', 'Weight output sensor [g]'),
                          ('gf_mt', 'Weight hanging sensor [g]')]:
        plt.figure(ind)
        ind += 1
        plt.ylabel(label)
        plt.xlabel('scan')

        value = eval(getattr(data, prefix))
        plt.plot(value, '.', label='measured', linewidth=2)

        for (suffix, leg) in [('linear', ['-', 'Lin. Sm.', 5]),
                              ('triangle',['-', 'Tri. Sm.', 5]),
                              ('gaussian', ['-', 'Gaus. Sm.', 5]),
                              ('gaussian', ['-', 'Gaus. Sm. deg 10', 10]),
                              ('gaussian', ['-', 'Gaus. Sm. deg 15', 15]),
                              ('gaussian', ['-', 'Gaus. Sm. deg 25', 25]),
                              ('gaussian_rec_2', ['-', 'Gaus. Sm. 2',5]),
                             ]:
            print (suffix.split('_'))
            if len(suffix.split('_')) > 1 and suffix.split('_')[1] == 'rec':
                sm_method = getattr(functions, 'smoothing_'+ suffix.split('_')[0]+ '_rec')
                plt.plot(sm_method(value, rec=int(suffix.split('_')[2]), degree=leg[2]), 
                         leg[0], label=leg[1],linewidth=2)
            else:
                sm_method = getattr(functions, 'smoothing_' + suffix)
                plt.plot(sm_method(value, degree=leg[2]), leg[0], label=leg[1], linewidth=2)

        legend = True
        plt.legend(loc=2)
        for (name, val) in zip(optname, optval):
            if name == 'xmin':
                plt.xlim(xmin=float(val))
            if name == 'xmax':
                plt.xlim(xmax=float(val))
            if name == 'legend':
                lval = int(val)
                legend = False
                if lval == -1:
                    pass
                else:
                    plt.legend(loc=lval)
        if legend:
            plt.legend()
    try:
        plt.show(block=False)
    except:
        plt.ion()
        plt.show()

if __name__ == "__main__":
    main()
    try:
        raw_input('Press key to quit')
    except:
        input('Press key to quit')
