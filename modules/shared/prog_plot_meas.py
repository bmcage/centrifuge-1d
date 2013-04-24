from __future__ import division, print_function

import numpy as np
import matplotlib.pyplot as plt

# 
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
    if len(ARGS) < 2:
        print ("ERROR: give ini file with gf meas data")
        sys.exit(0)
    fname = ARGS[1]
    if not os.path.isfile(fname):
        print ("ERROR: %s is not a file on your PC" % fname)
        sys.exit(0)

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
    for prefix in ['gf_mo', 'gf_mt']:
        plt.figure(ind)
        ind += 1
        plt.ylabel('weight %s' % prefix[-2:] )
        plt.xlabel('scan')
        for (suffix, leg) in [('', ['.','measured']), 
                               ('_sm', [':', 'Lin. Sm.']), 
                              ('_smtri',['-.', 'Tri. Sm.']), 
                              ('_smgau', ['-', 'Gaus. Sm.'])]:
            try:
                value = getattr(data, prefix+suffix)
                value = eval(value)
                plt.plot(value, leg[0], label=leg[1])
            except:
                pass
        plt.legend()
    try:
        plt.show(block=False)
    except:
        plt.ion()
        plt.show()

if __name__ == "__main__":
    main()
    raw_input('Press key to quit')