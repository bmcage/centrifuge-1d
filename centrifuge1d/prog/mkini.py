#!/usr/bin/python
from __future__ import print_function
import io, shutil
import csv
import numpy as np
from os import makedirs, path, sep
from optparse import OptionParser
from ..config import ModulesManager
from ..const import (CSV_DIR, INI_DIR, FIGS_DIR, DEFAULTS_ININAME,
                     MEASUREMENTS_ININAME, PROTO_DIR, PLOTSTYLE_ININAME)
from ..shared import get_directories
from collections import namedtuple, defaultdict
from ..modules.shared.measurements import MEASUREMENTS_NAMES
import logging
LOG = logging.getLogger(".")

MODMAN = ModulesManager()

BASE_MODULE = MODMAN.find_module('base', submodule='csv_parse')

# dynamicall created options: options added to generated MEASUREMENTS_ININAME
# file, that are not specified by the csv file
# {opt_name: gen_data_function} with 'gen_data_function(exp_type, data, modman)'

def parse_input():

    usage_str = (
         '\n\t %prog csv_ID'
         '\n\nwhere'
         '\n   csv_ID:'
         '\n       ID of the csv file to be processed. It will search for the '
         '.csv files in the CSV directory: \n       "' + CSV_DIR + '"'
         '\n       Inside has to be either file named csv_ID.csv (with .csv '
         'extension) or directory named \n       csv_ID containing file '
         'csv_ID.csv.'
         '\n\n       Example: '
         '\n           To process file named "experiment.csv" in the CSV '
         'directory (or in a subdirectory '
         '\n           named "experiment" inside the CSV directory), the '
         'command is:'
        '\n\t\t./mkini experiment')

    optparser = OptionParser(usage=usage_str)
    optparser.add_option('-l', '--list', dest='list', action="store_true",
                         default=False,
                         help="Lists all available CSV IDs")

    (options, args) = optparser.parse_args()

    len_args = len(args)

    if len_args == 0:
        if options.list:
            from os import listdir

            print('\nAvailable CSV IDs:\n--------------------')
            for csv_id in sorted(listdir(CSV_DIR)):
                if (len(csv_id) > 4) and (csv_id[-4:] == '.csv'):
                    print(csv_id[:-4])
                else:
                    print(csv_id)
            print()
        else:
            optparser.print_help()
        exit(0)

    csv_path = CSV_DIR + sep + args[0]

    if path.isdir(csv_path):
        csv_path        += sep
        csv_filebasename = args[0]
    else:
        csv_path         = CSV_DIR + sep
        csv_filebasename = args[0]

    if path.exists(csv_path + csv_filebasename):
        options.csv_filebasename = csv_filebasename[:-4] # truncate the extension
    elif path.exists(csv_path + csv_filebasename + '.csv'):
        options.csv_filebasename = csv_filebasename
    else:
        print('CSV file does not exist:{}'.format(csv_filebasename))
        exit(1)

    options.csv_path = csv_path

    if not path.exists(INI_DIR):
        makedirs(INI_DIR)

    return options

def get_module_function(function_name, module, base_module=BASE_MODULE):
    if hasattr(module, function_name):
        function = getattr(module, function_name)
    else:
        function = getattr(base_module, function_name)

    return function

def compress_seq(seq):
    if not type(seq) in [list, tuple]: return seq

    if not seq:
        return seq

    ref = seq[0]

    # all(seq[:] == ref)
    if not [value for value in seq if not (value == ref)]:
        return ref

    # all are == '' except one
    nonempty = [value for value in seq if not (value == '')]
    items_count = len(nonempty)
    if items_count == 1:
        return nonempty[0]
    elif items_count == len(seq):
        return seq
    else:
        raise ValueError('Either all values in sequence must be nonempty, '
                         'or exactly one:\n', seq)

############################################################
#          Storage class for data read from CSV file       #
############################################################
class ExperimentData():
    def __init__(self, data_type, initial_data={}):
        self._data_type = data_type
        self._data = defaultdict(list, initial_data)

    def get_value(self, key, not_found=None):
        if not key in self._data:
            return not_found
        else:
            return self._data[key]

    def set_value(self, key, value):
        self._data[key] = value

    def append_row(self, row):
        data = self._data

        for (key, value) in vars(row).items():
            data[key].append(value)

    def iterate(self, fn, user_data=None):
        if user_data is None:
            for (key, value) in self._data.items(): fn(key, value)
        else:
            for (key, value) in self._data.items(): fn(key, value, user_data)

    def delete(self, *keys):
        data = self._data
        for key in keys: del data[key]

    def on_modify(self, fn):
        """
        Set the value returned by function fn as new value. If new_value is
        'None', then the key is deleted. Type of fn: fn(key, value)->new_value
        """

        delete_items = []

        experiment_data = self._data

        for (key, value) in experiment_data.items():
            new_value = fn(key, value)
            if new_value is None: delete_items.append(key)
            else: experiment_data[key] = new_value

        self.delete(*delete_items)

    def save(self, save_dir, filename=DEFAULTS_ININAME,
             section_name='experiment', newlines=True):

        def _princ(key, value, stream):
            newlinecounter = 1
            addnewline = False
            if type(value) in (list, tuple, np.ndarray):
                if np.alen(value) == 0:
                    value_str = value
                else:
                    value_str = '[' + str(value[0])
                    for item in value[1:]:
                        if not addnewline:
                            value_str += ', ' + str(item)
                        else:
                            value_str += ',\n' + str(item)
                            addnewline = False
                        if newlines:
                            newlinecounter += 1
                            if newlinecounter > 9:
                                addnewline = True
                                newlinecounter = 0
                    value_str += ']'
            else:
                value_str = value

            stream.write('{:8} = {}\n'.format(key, value_str))

        if callable(save_dir):
            savedir = save_dir(self)
        else:
            savedir = save_dir

        if not path.exists(savedir):
            makedirs(savedir)

        try:
            fout = open(savedir + filename, mode='w', encoding='utf-8')
        except:
            fout = open(savedir + filename, mode='w')
        fout.write('[{}]\n'.format(section_name))

        self.iterate(_princ, user_data=fout)

        fout.close()

    def echo(self, short=False):
        print()
        if short:
            for (name, value) in sorted(self._data.items()):
                print('%-12s = %s' % ('key', name))
        else:
            for (name, value) in sorted(self._data.items()):
                print('%-12s = %s' % (name, value))

    def get_datatype(self):
        return self._data_type

class MultipleExperimentData():
    def __init__(self, data_type):
        self._data_type = data_type
        self._data = {}

    def pick(self):
        leaf = self._data

        while not isinstance(leaf, ExperimentData):
            leaf = leaf.keys()[0]

        return leaf

    def _get_leaf(self, ids, create_p=False, not_found=None):
        data          = self._data
        parental_data = None

        for identifier in ids:
            if not identifier in data:
                if create_p:
                    data[identifier] = {}
                else:
                    return not_found
            parental_data = data
            data = data[identifier]

        if data == {}:
            data = ExperimentData(self._data_type)
            parental_data[ids[-1]] = data

        return data

    def get_value(self, key, not_found=None, *ids):
        leaf = self._get_leaf(ids, not_found=EOFError)
        if leaf == EOFError:
            return not_found
        else:
            return leaf.get_value(key)

    def set_value(self, key, value, *ids):
        leaf = self._get_leaf(ids, create_p=True)
        leaf.set_value(key, value)

    def append_row(self, row, *ids):
        leaf = self._get_leaf(ids, create_p=True)
        leaf.append_row(row)

    def iterate(self, fn, user_data=None, level_prehook=None,
                level_posthook=None, levels_names=None):

        def _traverse(struct, level):
            if isinstance(struct, ExperimentData):
                if user_data is None:
                    fn(struct)
                else:
                    fn(struct, user_data)
            else:
                for (key, value) in struct.items():
                    if not level_prehook is None:
                        if levels_names:
                            if not user_data is None:
                                level_prehook(levels_names[level], key,
                                              user_data)
                            else:
                                level_prehook(levels_names[level], key)

                        else:
                            level_prehook(key, user_data)

                    _traverse(value, level+1)

                    if not level_posthook is None:
                        if levels_names:
                            if not user_data is None:
                                level_posthook(levels_names[level], key,
                                               user_data)
                            else:
                                level_posthook(levels_names[level], key)

                        else:
                            level_posthook(key, user_data)

        _traverse(self._data, 0)

    def delete(self, *keys):
        def _delete_on_leafs(experiment):
            experiment.delete(*keys)

        self.iterate(_delete_on_leafs)

    def on_modify(self, fn):
        def _modify_leafs(experiment):
            experiment.on_modify(fn)

        self.iterate(_modify_leafs)

    def save(self, save_dir, filename=DEFAULTS_ININAME,
             section_name='experiment'):

        def _save_configuration(experiment):
            experiment.save(save_dir, filename, section_name)

        self.iterate(_save_configuration, verbose=False)

    def get_datatype(self):
        return self._data_type

############################################################
#                      Read CSV file(s)                    #
############################################################
def missing_fields_in_csv(csv_fields, fields):

    missing = [field for field in fields if not field in csv_fields]

    return missing

def quote_string(string):
    tmp_str = string.strip()
    if len(tmp_str) > 1:
        if ((tmp_str[0] == tmp_str[-1] == '"')
            or (tmp_str[0] == tmp_str[-1] == "'")):
            pass
        else:
            tmp_str = "'" + tmp_str + "'"
    else:
        tmp_str = "'" + tmp_str + "'"

    return tmp_str

def quote_csv_file_string(key, value):

    if (key == 'csv_file') and (not value is None):
        new_value = []

        for csv_file_name in value:
            new_value.append(quote_string(csv_file_name))

        return new_value
    else:
        return value

def CSV2data(csv_dir, csv_filebasename):
    print('\t\t\tReading CSV file', csv_filebasename, '... ', end="")

    csv_filename = csv_dir + csv_filebasename + '.csv'
    # open CSV file
    try:
        f_csv = io.open(csv_filename, encoding='utf-8')
        csv_data = csv.reader(f_csv)
        csv_fields = next(csv_data)
    except UnicodeDecodeError:
        f_csv.close()

        try:
            f_csv = io.open(csv_filename, encoding='utf-16')
            csv_data = csv.reader(f_csv, delimiter=';')
            csv_fields = next(csv_data)
        except Exception:
            f_csv.close()

            print('\n\tCould not read CSV file:', csv_filebasename + '.csv',
                  'Unknown type of CSV file. Exiting...')
            exit(1)

    # determine CSV type
    if ((len(csv_fields) >= 2) and (csv_fields[0] == 'Name:')
        and (len(csv_fields[1]) >= 10) and (csv_fields[1][:10] == 'Data Instr')):

        try:
            #csv_fields = next(csv_data)
            while csv_fields[0] != 'Scan':
                csv_fields = next(csv_data) # skip useless header
        except Exception as msg:
            f_csv.close()

            print('Error during file reading occured. Message', msg, 'Exiting...')
            exit(1)

        # replace original column names with something more readable:
        # csv_fields = ('scan', 'time', 'elapsed', 's101', 's102', 's103', ...)
        dummy_nr = 1;
        for i in range(0, len(csv_fields)):
            if csv_fields[i] == '':
                csv_fields[i] = 'dummy' + str(dummy_nr)
                dummy_nr += 1
                continue

            csv_fields[i] = csv_fields[i].lower()
            fields = csv_fields[i].split(' ')
             # test for sensor measurements e.g. "101 (VDC)", "112 (VDC)"
            if (len(fields) > 1) and (fields[0][0] == "1"):
                csv_fields[i] = 's' + fields[0]

        csv_type = 'scanned_data'
        data = ExperimentData(csv_type)
    elif ((len(csv_fields) > 10) and (csv_fields[0] == 'nr')
        and (csv_fields[5] == 'ADC0')):

        csv_type = 'scanned_data_2017'
        data = ExperimentData(csv_type)

    else:
        csv_type = 'regular_ini'
        required_fields = ('exp_type', 'exp_id', 'exp_no')
        missing = missing_fields_in_csv(csv_fields, required_fields)
        if missing:
            if len(missing) == len(required_fields):
                print("Could not read csv table header. Possibly corrupted "
                      "csv file?")
            print("CSV file: missing identifiers: '%s' "% missing,
                  "\nCSV error: cannot process CSV file, aborting...")
            exit(1)

        extended_fields = ('scan', 'csv_file', 'omega')
        missing = missing_fields_in_csv(csv_fields, extended_fields)
        if not missing:
            csv_type = 'scan_ini'
        elif not 'scan' in missing:
            print("ERROR: 'scan' value is provided in CSV file, but needed "
                  "identifiers are missing: ",  missing)
            print('CSV error: cannot process CSV file, aborting...')
            exit(1)

        identifiers = {}

        data = MultipleExperimentData(csv_type)

    csv_row = namedtuple('CSVRow', csv_fields)

    try:
        for row in map(csv_row._make, csv_data):
            if csv_type in ('regular_ini', 'scan_ini'):
                module = MODMAN.find_module(row.exp_type, submodule='csv_parse')
                skip   = get_module_function('skip', module)

                if skip(row): continue

                for id_name in ('exp_type', 'exp_id', 'exp_no'):
                    value = getattr(row, id_name)
                    if value:
                        identifiers[id_name] = value # update if value present

                (exp_id, exp_no) = (identifiers['exp_id'], identifiers['exp_no'])
                data.append_row(row, exp_id, exp_no)
            else:
                data.append_row(row)
    except Exception as msg:
        f_csv.close()

        print('Error during file reading occured, msg=\n', msg, ' Exiting...')
        exit(1)

    if csv_type == 'scanned_data':
        for name in csv_fields:
            if name in ('s102', 's103'):
                # force to use '.' as a delimiter for floats (and not ',')
                dval = []
                for value in data.get_value(name):
                    dval.append(value.replace(',', '.'))

                data.set_value(name, dval)
            else:
                # remove useless columns from data structure
                data.delete(name)

    elif csv_type == 'scanned_data_2017':
        for name in csv_fields:
            if name in ('epochtime', 'kg_bottom', 'kg_hanging', 'RPM'):
                pass
            else:
                # remove useless columns from data structure
                data.delete(name)
    elif csv_type == 'scan_ini':
        # if needed, add quotes to values in csv_file field
        data.on_modify(quote_csv_file_string)

    return data

############################################################
#          Transform data and save to .ini file            #
############################################################
ITERABLE_LISTS = {}

def multiply_list(lst, coef):
    return [str(float(value) * coef) for value in lst]


def subtract_list(lst, term):
    return [str(float(value) - term) for value in lst]

def get_iterables_list(exp_type, modman):

    def extend_iterables_list(module):
        if hasattr(module, 'OPTIONS_ITERABLE_LISTS'):
            iterables_list.extend(module.OPTIONS_ITERABLE_LISTS)
        return True

    if exp_type in ITERABLE_LISTS:
        iterables_list = ITERABLE_LISTS[exp_type]
    else:
        iterables_list = []

        modman.traverse_ancestors(exp_type, extend_iterables_list,
                                  submodule='options')

        ITERABLE_LISTS[exp_type] = iterables_list

    return iterables_list

def find_first_idx(seq, start_index, treshold = 0.02, step = +1):
    if step > 0:
        last_index = len(seq) - 1
        test_index = lambda index: index <= last_index
    else:
        last_index = 0
        test_index = lambda index: index >= last_index


    value0 = float(seq[start_index])
    idx = start_index + step

    while test_index(idx):
        value1 = float(seq[idx])
        if (abs(value1 - value0) < treshold):
            break
        else:
            value0 = value1

        idx += step

    return idx

def determine_duration_and_omega(measurement_data, omegas_ini, scans, scan_span):

    c_omega = compress_seq(omegas_ini)

    if np.isscalar(c_omega): # omega is constant
        omega_cen = c_omega

        total_duration = scan_span * (scans[-1] - scans[0] + 1)
        include_acceleration = (not c_omega == 0.0)
        a_duration = include_acceleration * total_duration
        d_duration = 0.0
        g_duration = (not include_acceleration) * total_duration

    else: # omega is variable (we have (a+d)+g in general)
        omega_cen  = []
        a_duration = []
        d_duration = []
        g_duration = []

        omegas = np.asarray(omegas_ini, dtype=float)

        scans_nr = len(scans)

        omega_prev = omegas[0]
        scan_prev  = scans[0]

        if omegas[0] == 0.0:
            a_duration.append(0.0)
            d_duration.append(0.0)
            omega_cen.append(0.0)

        for i in range(1, scans_nr):
            omega = omegas[i]
            scan = scans[i]

            if (omega == omega_prev) and (i+1 < scans_nr):
                continue

            if omega_prev == 0.0:    # acceleration phase
                g_duration.append(scan_span * (scan - scan_prev))

            elif omega == 0.0:       # g-phase
                omega_cen.append(omega_prev)
                a_duration.append(scan_span * (scan - scan_prev))
                d_scan = find_first_idx(measurement_data, scan,
                                        step=-1, treshold=0.02)
                d_duration.append(scan_span * (scan - d_scan))

            elif omega_prev <= omega: # (continuation of) acceleration phase
                # omega_prev == omega <=> i+1 == scans_nr
                omega_cen.append(omega_prev)
                a_duration.append(scan_span * (scan - scan_prev))
                d_duration.append(0.0)
                g_duration.append(0.0)

            else: # 0 < omega_next < omega
                    raise NotImplementedError('Decelerating to omega>0 is not '
                                              'supported.')
            omega_prev = omega
            scan_prev  = scan

        if (omega == 0.0) and (omegas[-2] > 0):
            g_duration.append(0.0)
        elif (omega > 0.0) and (not omegas[-2] == omega):
            print('\nWARNING: Last omega must be of value either 0.0 '
                  'or it''s previous value. Assuming the latter.\n')

        if compress_seq(a_duration) == 0.0: a_duration = 0.0
        if compress_seq(d_duration) == 0.0: d_duration = 0.0
        if compress_seq(g_duration) == 0.0: g_duration = 0.0

    return (omega_cen, a_duration, d_duration, g_duration)

CSV2INI_NAMES = {'s102': 'gf_mt', 's103': 'gf_mo', 'kg_bottom': 'gf_mo',
                 'kg_hanging': 'gf_mt', }

def process_scanned_data(experiment, experiment_info):
    csv_files = experiment.get_value('csv_file')
    if not csv_files or csv_files in ["''", '""']:
        return

    data_dir  = get_directories('ini', 'data', experiment_info)

    if type(csv_files) == str:
        csv_scanned = csv_files.strip()
        if ((csv_scanned[0] == csv_scanned[-1] == '"')
            or (csv_scanned[0] == csv_scanned[-1] == "'")):
            csv_scanned = csv_scanned[1:-1]
        scanned_data = CSV2data(experiment_info['csv_path'], csv_scanned)
    else:
        raise ValueError('CSV file describing scanned data can contain '
                         'only one filename: ', csv_files)
    omegas = experiment.get_value('omega')
    scans  = np.asarray(experiment.get_value('scan'), dtype=int)
    scan_span = float(experiment.get_value('scan_span', not_found=1.0))

    if scanned_data.get_datatype() == 'scanned_data':

        scanned_inflow_weights  = scanned_data.get_value('s102')
        if scans[-1] == '':
            scans[-1] = len(scanned_inflow_weights) - 1
        (scan_first, scan_last) = (scans[0], scans[-1])

        if not omegas:
            print('Rotational speed ''omega'' not specified: ', omegas,
                  '\nExiting...')
            exit(1)
        (omega_cen, a_duration, d_duration, g_duration) = \
          determine_duration_and_omega(scanned_inflow_weights, omegas, scans,
                                       scan_span)
        include_acceleration = (not d_duration == 0.0)

        experiment.set_value('omega', omega_cen)
        experiment.set_value('include_acceleration', include_acceleration)
        experiment.set_value('duration', a_duration)
        experiment.set_value('deceleration_duration', d_duration)
        experiment.set_value('fh_duration', g_duration)

        # Truncate unused values and store xvalues
        times = scan_span * np.arange(1, scan_last-scan_first+1)
        for (csv_name, ini_name) in CSV2INI_NAMES.items():

            # truncate
            scanned_force_weights = scanned_data.get_value(csv_name)
            if scanned_force_weights:
                # store xvalue
                scanned_data.set_value(ini_name + '_xvalues', times)
                scanned_force_weights = \
                  multiply_list(scanned_force_weights, 1000.) # kg -> g
                scanned_data.set_value(ini_name,
                                       scanned_force_weights[scan_first: scan_last])
                scanned_data.delete(csv_name)

    elif scanned_data.get_datatype() == 'scanned_data_2017':

        scanned_inflow_weights  = scanned_data.get_value('kg_hanging')
        if scans[-1] == '':
            scans[-1] = len(scanned_inflow_weights) - 1
        (scan_first, scan_last) = (scans[0], scans[-1])

        # deterime time of experiments
        times = scanned_data.get_value('epochtime')
        times = times[scan_first: scan_last]
        starttime = float(times[0])
        times = subtract_list(times, starttime)
        scanned_data.delete('epochtime')

        # Truncate unused values and store xvalues
        for (csv_name, ini_name) in CSV2INI_NAMES.items():
            # truncate
            scanned_force_weights = scanned_data.get_value(csv_name)
            if scanned_force_weights:
                # store xvalue
                scanned_data.set_value(ini_name + '_xvalues', times)
                scanned_force_weights = \
                  multiply_list(scanned_force_weights, 1000.) # kg -> g
                scanned_data.set_value(ini_name,
                                       scanned_force_weights[scan_first: scan_last])
                scanned_data.delete(csv_name)
#
#        #remove unneeded data that was scanned
#        for scanheader in ['nr', 'date', 'time', 'gain', 'ADC0', 'ADC1', 'ADC2',
#                           'ADC3', 'V0', 'V1', 'V2', 'V3', 'frqHz',
#                           'frqRPM' ]:
#            scanned_data.delete(scanheader)

        if not omegas:
            print('Rotational speed ''omega'' not specified: ', omegas,
                  '\nExiting...')
            exit(1)

        c_omega = compress_seq(omegas)

        if np.isscalar(c_omega): # omega is constant
            omega_cen = [omegas]

        else: # omega is variable
            omega_cen = np.asarray(omegas, dtype=float)
            if '-1' in omega_cen:
                print('Rotational speed ''omega'' -1 present, but not -1 for sublines. All -1 or None! ', omegas,
                     '\nExiting...')
            exit(1)

        if omega_cen[0] != '-1':
            # use generated omega data
            (omega_cen, a_duration, d_duration, g_duration) = \
              determine_duration_and_omega(scanned_inflow_weights, omegas, scans,
                                           scan_span)
            include_acceleration = (not d_duration == 0.0)

            experiment.set_value('omega', omega_cen)
            experiment.set_value('include_acceleration', include_acceleration)
            experiment.set_value('duration', a_duration)
            experiment.set_value('deceleration_duration', d_duration)
            experiment.set_value('fh_duration', g_duration)
        else:
            #use omega as measured by rpm sensor
            scanned_omega = scanned_data.get_value('RPM')
            if scanned_omega:
                # store xvalue
                scanned_data.set_value('radps' + '_xvalues', times)
                scanned_omega = \
                  multiply_list(scanned_omega, 2*np.pi / 60.) # rpm -> rad/s
                scanned_data.set_value('radps',
                                       scanned_omega[scan_first: scan_last])
                scanned_data.delete('RPM')
            else:
                print (scanned_data.echo(short=True))
                raise ValueError('RPM data not present in scanned data, but required.')

    else:
        raise ValueError('CSV file of type', scanned_data.get_datatype(),
                         'which is an unknown data type')

    # Remove unused data
    experiment.delete('scan')

    filename = csv_scanned + '.ini'
    scanned_data.save(data_dir, filename)

PROCESS_DATA = {'regular_ini': lambda exp, info: True,
                'scan_ini': process_scanned_data}

def data2ini(data, experiment_info):
    def _filter_data(key, value):
        cvalue = compress_seq(value)

        # durations need to be handled specially
        if key in ('duration', 'fh_duration', 'deceleration_duration'):
            if ((not type(cvalue) in (list, tuple))
                and float(cvalue) == 0.0):

                return 0.0
            else: return value

        if cvalue in [None, '']: return None
        else: return cvalue

    def _pre_process_data(lname, lvalue, experiment_info):
        if lname == 'exp_id':
            print("\tProcessing experiment(s) with ID: ", repr(lvalue))
        else:
            print("\t\tProcessing experiment number: ", lvalue)

    def _post_process_data(lname, lvalue, experiment_info):

        if lname == 'exp_id':
            print("\tCreating default inifile for ID ",
                  repr(experiment_info['exp_id']), end = '')

            exp_base_path = get_directories('ini', 'exp_base', experiment_info)

            if path.exists(exp_base_path + DEFAULTS_ININAME):
                print('\n\t\tFile {0} exists.\t\t\t\t\tSkipping'
                      .format(exp_base_path + DEFAULTS_ININAME))

            else:
                default_configuration = \
                  ExperimentData('defaults',
                                 {'exp_type': repr(experiment_info['exp_type'])})

                default_configuration.save(exp_base_path,
                                           filename=DEFAULTS_ININAME,
                                           section_name='experiment')

            if not path.exists(exp_base_path + PLOTSTYLE_ININAME):
                LOG.warn ("copying " + PROTO_DIR + PLOTSTYLE_ININAME + '.proto' + ' to ' + exp_base_path + PLOTSTYLE_ININAME)
                shutil.copy(PROTO_DIR + PLOTSTYLE_ININAME + '.proto',
                            exp_base_path + PLOTSTYLE_ININAME)

            # copy the original csv files (in case they are overwritten in the
            # future so that we know what the data was generated from)
            print('\nSaving the original csv file(s)...')
            csv_filename = experiment_info['csv_filebasename'] + '.csv'
            shutil.copy(experiment_info['csv_path'] + csv_filename,
                        exp_base_path + 'orig.' + csv_filename)

            experiment_info['exp_type'] = None # reset exp_id for _process_data

            print()

    def _process_data(experiment, experiment_info):
        # Process single configuration file
        experiment.on_modify(_filter_data) # 'None' value returned from
                                           # _filter_data() is deleted

        exp_type = experiment.get_value('exp_type')
        if type(exp_type) in (list, tuple):
            print("Experiment type 'exp_type' must be the "
                  "same for given experiment. Exiting...")
            exit(1)

        if experiment_info['exp_type'] is None:
            experiment_info['exp_type'] = exp_type
        elif experiment_info['exp_type'] != exp_type:
            print("Experiment type 'exp_type' must be the same "
                      "for all experiments with given 'exp_id'.")
            exit(1)

        experiment_info['exp_id'] = experiment.get_value('exp_id')
        experiment_info['exp_no'] = experiment.get_value('exp_no')

        (data_dir, masks_dir)  = \
          get_directories('ini', ['data', 'masks'], experiment_info)

        if not path.exists(masks_dir):
            makedirs(masks_dir) # create empty masks directory

        # support custom modification per-data_type
        PROCESS_DATA[experiment.get_datatype()](experiment, experiment_info)

        experiment.delete('exp_id', 'exp_no', 'exp_type')

        experiment.save(data_dir, filename=MEASUREMENTS_ININAME,
                        section_name='experiment-data')

        print('\t\t\tDone')

    print('\nCreating data inifiles...')
    experiment_info.update(dict.fromkeys(('exp_id', 'exp_no', 'mask',
                                          'exp_type')))
    data.iterate(_process_data, user_data=experiment_info,
                 levels_names=('exp_id', 'exp_no'),
                 level_prehook=_pre_process_data,
                 level_posthook=_post_process_data)
    print('All done.')

def main(csv_path=None, csv_filebasename=None, dataout_path=None):
    """
    Make ini files based on  a file.csv in csv_path. Here file would be
    csv_filebasename.
    If one of the arguments is not given, argv is parsed to determine them.
    """
    if csv_path is None or csv_filebasename is None:
        options = parse_input()
        info = {'csv_path': options.csv_path,
               'csv_filebasename': options.csv_filebasename}
        if csv_path:
            info['csv_path'] = csv_path
        if csv_filebasename:
            info['csv_filebasename'] = csv_filebasename
    else:
        info = {'csv_path': csv_path + sep,
               'csv_filebasename': csv_filebasename}

    if dataout_path:
        info['ini_dir'] = dataout_path + sep + 'datafiles'
        info['figs_dir'] = dataout_path + sep + 'datafiles'
    else:
        info['ini_dir'] = INI_DIR
        info['figs_dir'] = FIGS_DIR

    data = CSV2data(info['csv_path'], info['csv_filebasename'])
    data2ini(data, info)

if __name__ == "__main__":
    main()
