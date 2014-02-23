from __future__ import print_function, division

import numpy as np, matplotlib, pickle
import sys
from ...const import PLOTSTYLE_ININAME, DUMP_DATA_VERSION, DUMP_DATA_FILENAME
from os import makedirs, path
from ...shared import get_directories, parse_value, yn_prompt, filter_indices
from ...config import ModulesManager, load_model
from collections import OrderedDict
from .functions import has_data, compare_data
try:
    import ConfigParser as configparser
except:
    import configparser

################################################################################
#                   Define graphs, axes labels and units                       #
################################################################################

DATA_UNITS = {'length': ('mm', 'cm', 'm'),
              'time': ('s', 'min', 'h'),
              'pressure': ('Pa', 'kPa'),
              'weight': ('g', 'kg'),
              'force_kgp': ('gf', 'kgf'),
              'velocity': ('cm/s', 'm/s'),
              'rotational_speed': ('rpm', 'rps', 'rad/s', 'rad/min'),
              'none': ('', )}
DEFAULT_UNITS = {'length': 'cm', 'time': 'min', 'pressure': 'Pa',
                 'force_kgp': 'gf', 'weight': 'g', 'velocity': 'cm/s',
                 'rotational_speed': 'rpm', 'none': ''}

# Linestyle is searched in order:
#   1. Exists 'lineid' and there for 'fig_id' specific option
#   2. If not, search '_default_' for 'fig_id'
#   3. If not, use option in the '_base_' of the '_default_'
# Linestyle options:
#    Assigned only to '_base_' of line_id:
#        'xdata', 'ydata', 'label', 'legend_data', 'order',
#    Assigned anywhere:
#        'width', 'symbolsize', 'lineopt', 'ls'
LINESTYLES_DEFAULT = {'_default_': {'_base_': {'width': 1, 'symbolsize': None,
                                               'lineopt': '.', 'order': 999},
                                    'h':      {'lineopt': '-'},
                                    'u':      {'lineopt': '-'},
                                    'theta':  {'lineopt': '-'},
                                    'relsat': {'lineopt': '-'}},
                      'measured':  {'_base_': {'lineopt': 'x'},
                                    'theta':  {'lineopt': '.'},
                                    'relsat': {'lineopt': '.'}},
                      'original':  {'_base_': {'lineopt': '-'}}
                     }

dg_label_time = "Time [{}]"
dg_label_length = "Sample length $L$ [{}]"
dg_unit_time_length = ('time', 'length')

# Figures options:
#    'order', 'xlabel', 'ylabel',
#    'xscale', 'yscale', 'xunit', 'yunit', 'xmin', 'ymin', 'xmax', 'ymax',
#    'show', 'ls', 'show_legend', 'legend_title', 'legend_bbox', 'legend_loc'

FIG_OPTIONS_DEFAULTS = {'show': True, 'legend_loc': 4, 'order': 999,
                        'show_legend': True}


FIGURES_DEFAULTS = \
              {'h':  {'title':  'Hydraulic head',
                      'xlabel': dg_label_length,
                      'ylabel': "Piezometric head $h$ [{}]",
                      'xtype':  'length',
                      'ytype':  'length',
                      'legend_title': 'Time [min]',
                      'legend_bbox': (1.01, 1.),
                      'legend_loc': 2},
               'u':  {'title': 'Effective saturation',
                      'xlabel': dg_label_length,
                      'ylabel': "Effective saturation $S_e${}",
                      'xtype':  'length',
                      'ytype':  'none',
                      'legend_title': 'Time [min]',
                      'legend_bbox': (1.01, 1.),
                      'legend_loc': 2},
               'MO': {'title':  'Amount of expelled water',
                      'xlabel': dg_label_time,
                      'ylabel': "Expelled water [{}]",
                      'xtype':  'time',
                      'ytype':  'length'},
               'MI': {'title':  'Amount of water in the inflow chamber',
                      'xlabel': (dg_label_time, "Inflow water [{}]"),
                      'xtype':  'time',
                      'ytype':  'length'},
               'GC': {'title':  'Gravitational centre',
                      'xlabel': dg_label_time,
                      'ylabel': "Gravitational centre [{}]",
                      'xtype':  'time',
                      'ytype':  'length'},
               'RM': {'title':  'Rotational momentum',
                      'xlabel': dg_label_time,
                      'ylabel': "Rotational momentum [$kg.m/s$]{}",
                      'xtype':  'time',
                      'ytype':  'none'},
               's1': {'title':  'Interface s1',
                      'xlabel': dg_label_time,
                      'ylabel': "Interface s1 [{}]",
                      'xtype':  'time',
                      'ytype':  'length'},
               's2': {'title':  'Interface s2',
                      'xlabel': dg_label_time,
                      'ylabel': "Interface s2 [{}]",
                      'xtype':  'time',
                      'ytype':  'length'},
               'WM': {'title':  'Water mass balance',
                      'xlabel': (dg_label_time, "Water mass balance [{}]"),
                      'xtype':  'time',
                      'ytype':  'length'},
               'WM_in_tube':  {'title':       '"Water mass in tube',
                               'xlabel': dg_label_time,
                               'ylabel': "Water mass in tube [{}]",
                               'xtype':  'time',
                               'ytype':  'length'},
               'theta':   {'title':       'Water content $\\theta${}',
                           'xlabel': "Water content $\\theta${}",
                           'ylabel': "Pressure $p$ [{}]",
                           'xtype':  'none',
                           'ytype':  'pressure',
                           'yscale': 'log',
                           'legend_loc': 1},
               'relsat':  {'title':       'Effective saturation $S_e$',
                           'xlabel': "Effective saturation $S_e${}",
                           'ylabel': "Negative hydraulic head $h$ [{}]",
                           'xtype':  'none',
                           'ytype':  'length',
                           'yscale': 'log',
                           'legend_loc': 1},
               'K':       {'title':  'Water content $\\theta$',
                           'xlabel': "Water content $\\theta${}",
                           'ylabel': "Hydraulic conductivity $K(\\theta)$ [{}]",
                           'xtype':  'none',
                           'ytype':  'velocity',
                           'yscale': 'log'},
               'K_u':     {'title':  'Effective saturation',
                           'xlabel': "Effective saturation $S_e${}",
                           'ylabel': "Hydraulic conductivity $K(S_e)$ [{}]",
                           'xtype':  'none',
                           'ytype':  'velocity',
                           'yscale': 'log'},
               'gF_MO':   {'title':  'Force of expelled water',
                           'xlabel': dg_label_time,
                           'ylabel':  "Force of expelled water [{}]",
                           'xtype':  'time',
                           'ytype':  'force_kgp'},
               'gF_MT':   {'title':  'Force of water in tube',
                           'xlabel': dg_label_time,
                           'ylabel': "Force of water in tube [{}]",
                           'xtype':  'time',
                           'ytype':  'force_kgp'},
               'dgF_MO':  {'title':  'Force difference of expelled water',
                           'xlabel': dg_label_time,
                           'ylabel':  "Force difference of expelled water [{}]",
                           'xtype':  'time',
                           'ytype':  'force_kgp'},
               'dgF_MT':  {'title':  'Force difference of water in tube',
                           'xlabel': dg_label_time,
                           'ylabel':  "Force difference of water in tube [{}]",
                           'xtype':  'time',
                           'ytype':  'force_kgp'},
               'omega_rpm': {'title':  'Rotational speed',
                             'xlabel': dg_label_time,
                             'ylabel':  "Rotational speed [{}]",
                             'xtype':  'time',
                             'ytype':   'rotational_speed'}
                      }

FIGURES_PAIRS = (('h', 'u'), ('MI', 'MO'), ('GC', 'RM'), ('thK', 'uK'),
                 ('gF_MT', 'gF_MO'), ('dgF_MT', 'dgF_MO'), ('s1', 's2'),
                 ('theta', 'relsat'))

FIGURES_IDS = list(FIGURES_DEFAULTS.keys())

DISPLAY_OPTIONS = {'separate_figures': False, 'show_figures': True,
                   'save_figures': True, 'matplotlib_backend': None}

def set_default_units(figures_styles):
    for (fig_id, fig_style) in figures_styles.items():
        fig_style['xunit'] = DEFAULT_UNITS[fig_style['xtype']]
        fig_style['yunit'] = DEFAULT_UNITS[fig_style['ytype']]

    return figures_styles

def get_unit_coef(unit_base):
    """
      Return the coeficient for the new unit type to be used, so that
      the internal unit is converted to the new unit.
    """

    unit = unit_base.lower()
    # units used for computation are: cm, s, pa, gf and "no units"
    if unit in ['cm', 's', 'pa', 'g', 'gf', 'cm/s', 'rpm', '']: coef = 1.0
    elif unit == 'mm': coef = 10.
    elif unit in ('min'): coef = 1./60.
    elif unit == 'h': coef = 1./3600.
    elif unit in ('m', 'm/s'): coef = 0.01
    elif unit in ['kpa', 'kg', 'kgp']: coef = 0.001
    elif unit == 'rps': coef = 60.
    elif unit == 'rad/s': coef = np.pi/30.
    elif unit == 'rad/min': coef = 2*np.pi
    else:
        print('Unknown unit:', unit_base, '\nKnown units are only:', DATA_UNITS)
        exit(1)
    return coef

################################################################################
#                   Methods for handling default values                        #
################################################################################

def get_line_option(linestyles, line_id, name, fig_id=None, not_found=None):
    if fig_id is None:
        # we search only values in '_base_'
        if (line_id in linestyles) and (name in linestyles[line_id]['_base_']):
            value = linestyles[line_id]['_base_'][name]
        elif name in linestyles['_default_']['_base_']:
            value = linestyles['_default_']['_base_'][name]
        else:
            value = not_found

        return value

    # search for most specific to least specific
    for line_name in (line_id, '_default_'):
        if not line_name in linestyles:
            continue

        line_style = linestyles[line_name]

        if (fig_id in line_style) and (name in line_style[fig_id]):
            value = line_style[fig_id][name]
            return value

        elif name in line_style['_base_']:
            value = line_style['_base_'][name]
            return value

    return not_found

def get_figure_option(figstyles, fig_id, name, not_found = None):
    fig_style = figstyles[fig_id]

    if name in fig_style:
        value = fig_style[name]
    elif name in FIG_OPTIONS_DEFAULTS:
        value = FIG_OPTIONS_DEFAULTS[name]
    else:
        value = not_found

    return value

################################################################################
#                                Data storage                                  #
################################################################################

class DataStorage():
    """
      Object for holding computed, measured data and additional data.
    """

    def __init__(self, experiment_info):
        self._data = {'lines': {},
                      'experiment_info': experiment_info}
        self._modman = None

    def has_data(self, data_type):
        """ Return true if data 'data_type' is present in stored data. """

        return ((data_type in self._data)
                and (not self._data[data_type] is None))

    def store(self, name, value):
        """ Store the value under supplied name. """
        self._data[name] = value

    def get(self, key, not_found=None):
        """ Get the value corresponding to supplied name. """
        if key in self._data:
            return self._data[key]
        else:
            return not_found

    def store_computation(self, model, measurements, ID='computed'):
        """ Store computed values with ID """
        if self._modman is None:
            self._modman = ModulesManager()

        solver_module = self._modman.find_module(model.exp_type, submodule='run')
        if not hasattr(solver_module, 'solve_direct'):
            print("WARNING: Module for type '" + model.exp_type + "' does not "
                  "specify a 'solve_direct' variable. Function 'solve()' will "
                  "be used instead.")
            flag = solver_module.solve(model, measurements)
        else:
            flag = solver_module.solve_direct(model, measurements)

        if not flag:
            print("Computation of reference with ID '" + ID + "' was not "
                  "successfull. Data will not be saved.")
        else:
            data = {}

            # Store all computed data
            for (name, xvalue, yvalue) in measurements.iterate_calc_measurements():
                if xvalue is None or yvalue is None:
                    continue
                # make a local copy as array may be overwritten
                xvalue = xvalue.copy()
                yvalue = yvalue.copy()

                if name in ('h', 'u'):
                    t = measurements.get_times()
                    data[name] = (xvalue.transpose(), yvalue.transpose(), t)
                else:
                    data[name] = (xvalue, yvalue)

            # Store extra data
            # a) Retention curve based on theta
            if hasattr(model, 'SC'):

                SC = model.SC

                if hasattr(model, 'theta_s'): theta_s = model.theta_s
                else: theta_s = model.porosity

                if hasattr(model, 'theta_r'): theta_r = model.theta_r
                else: theta_r = 0.0

                (p, theta) = SC.retention_curve(theta_s, model.density, model.g,
                                                theta_r=theta_r)
                (rc_h, rc_u) = SC.retention_curve_u(g=model.g,
                                                    rho=model.density)

                if hasattr(model, 'p') or hasattr(model, 'pc'):
                    if hasattr(model, 'p'):
                        p_meas = np.asarray(model.p)
                    else:
                        p_meas = np.asarray(model.pc)

                    (p_user, theta_user) = \
                      SC.retention_curve(theta_s, model.density, model.g,
                                         theta_r=theta_r, p=p_meas)
                    (rc_h_user, rc_u_user) = \
                      SC.retention_curve_u(g=model.g, rho=model.density,
                                           p=p_meas)

                    data['theta']  = (p, theta, p_user, theta_user)
                    data['relsat'] = (-rc_h, rc_u, -rc_h_user, rc_u_user)
                else:
                    data['theta']  = (p, theta)
                    data['relsat'] = (-rc_h, rc_u)

                if hasattr(model, 'ks'):
                    try:
                        thK = SC.conductivity_curve(model.ks, theta_s,
                                              theta_r=theta_r, g=model.g,
                                              rho=model.density, rtype='theta')

                        uK = SC.conductivity_curve(model.ks, theta_s,
                                              theta_r=theta_r, g=model.g,
                                              rho=model.density, rtype='u')
                        data['thK'] = thK
                        data['uK']  = uK
                    except:
                        import traceback
                        print (traceback.format_exc())

            self._data['lines'][ID] = data
            self.store('experiment_info', model.experiment_info)

        return flag

    def store_references(self, user_references, model=None):
        """ Store computations corresponding to referenced parameters. """

        stored_references = self.get('references', not_found={})

        if type(user_references) in (list, tuple):
            print("DEPRECATION ERROR: format of 'params_ref' has been changed."
                  "Please update to new format.")
            return False

        if user_references == stored_references:
            return False   # nothing needs to be re-computed

        if model is None:
            model = load_model(self.get('experiment_info'), validate=True)
            #set original parameters as the found inverse param
            model.set_parameters(self.get('inv_params', {}))

        for ref_id in list(stored_references.keys()): # remove nonexisting refs
            if not ref_id in user_references:
                del stored_references[ref_id]

        iterable_params =  model._iterable_parameters
        for (ref_id, ref_params) in user_references.items():
            if ((ref_id in stored_references) # stored value is the same
                and (stored_references[ref_id] == ref_params)):
                continue

            ref_params = ref_params.copy() # work with backup
            stored_references[ref_id] = ref_params

            iters = [val for val in ref_params if val in iterable_params]
            if iters:
                print('Referencing model cannot set iterable '
                      'parameters of original model:', iters,
                      '\nSkipping...')
                continue

            backup_params = model.get_parameters(ref_params) # backup
            if hasattr(model, 'SC'):
                backup_SC = model.SC
                backup_typeSC = model.SC.typeSC()
                thekey = [key.lower() for key in ref_params.keys()]
                if not 'sc_type' in thekey:
                    print ('Referencing model does not contain "SC_type", cannot '
                           'set parameters of the saturation curve.'
                           '\nSkipping...')

            model.set_parameters(ref_params)

            flag = self.store_computation(model, model.measurements, ID=ref_id)

            if not flag:
                print('Reference parameters: ', ref_params)

             # restore original SC
            if hasattr(model, 'SC'):
                if backup_typeSC !=  model.SC.typeSC():
                    model.SC = backup_SC
            model.set_parameters(backup_params)

        self.store('references', stored_references)

        return True

    def store_measurements(self, measurements, model=None):
        """ Store measured (supplied) data. """

        for mtype in ('original', 'measured'):
            m = {}

            untransformed = (mtype == 'original')
            for (name, xvalue, yvalue) \
                in measurements.iterate_meas_measurements(untransformed, model):

                m[name] = (xvalue, yvalue)

                self._data['lines'][mtype] = m

    def get_linedata(self, line_id, not_found=None):
        """ Get data for specified line. """

        data = self._data['lines']

        if line_id in data:
            return data[line_id]
        else:
            return not_found

    def save(self):
        """ Save all stored data to a file. """

        if not self._data:
            print('No data was stored. Nothing to be saved. Skipping saving...')
            return

        savedir = get_directories('figs', 'mask', self.get('experiment_info'))
        if not path.exists(savedir):
            makedirs(savedir)

        with open(savedir + DUMP_DATA_FILENAME, 'wb') as fout:
            pickle.dump(self._data, fout, DUMP_DATA_VERSION)

    def load(self):
        """ Load data stored in a file. """

        pathdir = get_directories('figs', 'mask', self.get('experiment_info'))
        filename = pathdir + DUMP_DATA_FILENAME
        if not path.exists(filename):
            print('INFO: File with computation results does not exist:',
                  filename)
            return False

        with open(filename, 'rb') as fout:
            self._data = pickle.load(fout)
            # Old compatibility
            if 'ResultsData' in self._data:
                self._data = self._data['ResultsData']

        return True

################################################################################
#          Reading display configuration files and data displaying             #
################################################################################

def print_status(data, filename=None):
    """
      Compare measured vs. computed data.
      If filename is supplied, write it to the file.
    """

    computed     = data.get_linedata('computed')
    measurements = data.get_linedata('measured')

    if not computed: return

    if filename is None:
        stream = None
    else:
        stream = open(filename, 'w')

    if not measurements: return

    try:
        # compare measured vs. computed data
        for (key, m_data) in measurements.items():
            if m_data[1] is None: continue

            if key in computed:
                if key == 'theta':
                    c_value = computed[key][3]
                else:
                    c_value = computed[key][1]
                m_value = m_data[1]
            else:
                continue

            if c_value is not None:
                compare_data(key, c_value, m_value, stream)

        # display additional data
        cov = data.get('cov')
        if not cov is None:
            print('\nCov:\n', cov, file=stream)
            print('Deviation:\n', np.sqrt(np.diag(cov)), file=stream)

        params = data.get('inv_params')
        if not params is None:
            print('\nOptimal parameters found:', file=stream)
            for (name, value) in params.items():
                if name == 'ks':
                    print('  Ks [cm/s]: {: .8g}'.format(value), file=stream)
                else:
                    if np.iterable(value):
                        print('  {0:9}: {1!s}'.format(name, value), file=stream)
                    else:
                        print('  {:9}: {: .8g}'.format(name, value), file=stream)

    finally:
        if not filename is None: stream.close()

def get_filenames(experiment_info):
    """ Given experiment information find a list of all plotstyle inifiles. """

    (search_dirs, masks_dir, mask_dir) = \
        get_directories('figs', ['search', 'masks', 'mask'], experiment_info)

    search_dirs.append(masks_dir) # search also in masks_dir...
    if experiment_info['mask']:   # and if exists the mask directory, there too
        search_dirs.append(mask_dir)

    plotstyles_files = \
      [fname for fname in [cfgdir + PLOTSTYLE_ININAME for cfgdir in search_dirs]
       if path.exists(fname)]

    return plotstyles_files

def deep_dictupdate(d1, d2):
    """ Recursively update dictionary d1 with values of dictionary d2. """

    for (vkey, vvalue) in d2.items():
        if ((vkey in d1) and (type(d1[vkey]) == dict)
            and (type(vvalue) == dict)):
            deep_dictupdate(d1[vkey], vvalue)
        else:
            d1[vkey] = vvalue

    return d1

def read_plotstyles_file(filename):
    """ Reads the plotstyles configuration from files. """

    result = {}

    # read config files
    parser   = configparser.ConfigParser()
    try:
        parser.read(filename)
    except configparser.ParsingError as E:
        print('Plotstyles file ', filename,
              '\n could not be read because of the following error: ', E,
              '\nDefault values are used.')

        return result

    # Write data from parser to configuration
    for psection in parser.sections():
        # we ignore sections
        try:
            for option in parser.options(psection):
                raw_value = parser.get(psection, option).strip()

                result[option] = parse_value(raw_value, raise_error=True)
        except:
            if yn_prompt('WARNING: Errors during parsing occured. Do you wish '
                         'to continue with default values used? [Y/n]: '):
                result = {} # whatever was read so far is discared
                break
            else:
                exit(0)

    return result

def update_styles(styles, user_styles):
    """
      Update the styles from user plotstyles.ini file.
      Perform also any necessary conversion to the internal format.
    """

    # process 'lines' separately - just put missing values into
    # the default lines object
    user_lines_styles = user_styles['lines']
    lines_styles = styles['lines']

    for (line_id, data) in user_lines_styles.items():
        if not line_id in lines_styles:
            lines_styles[line_id] = {}

        line_style = lines_styles[line_id]
        if not '_base_' in line_style:
            line_style['_base_'] = {}

        for (name, value) in data.items():
            if name in FIGURES_IDS:
                # value is lineopt
                if not name in line_style:
                    line_style[name] = {}

                line_style[name]['lineopt'] = value
            else:
                # value is something specific to whole line, e.g. label
                line_style['_base_'][name] = value

    # Don't update with the remaining values
    del user_styles['lines']

    # Update all the rest
    styles = deep_dictupdate(styles, user_styles)

    # Display legend only on 'u' if 'u' and 'h' are on the same window
    styles['figures']['h']['show_legend'] = \
      styles['options']['separate_figures']

    # Order lines
    lines_ids = (['original', 'measured', 'computed']
                 + list(styles['params_ref'].keys()))
    lines_order = {line_id: get_line_option(lines_styles, line_id, 'order')
                   for line_id in lines_ids}
    styles['lines_order'] = \
      list(sorted(lines_order, key=lines_order.__getitem__))

    return styles

def order_figures(figures_styles, figs_ids):
    """
      Sort figures identified by 'figs_ids' according to their 'order'
      value. Pairs defined by FIGURES_PAIRS take precedence.
    """

    ordered_figures = []

    if not figs_ids:
        return ordered_figures

    figs_order = {fig_id: get_figure_option(figures_styles, fig_id, 'order')
                   for fig_id in figs_ids}

    for (id1, id2) in FIGURES_PAIRS:
        if (id1 in figs_order) and (id2 in figs_order):
            ordered_figures.append(id1)
            ordered_figures.append(id2)

            del figs_order[id1]
            del figs_order[id2]

    ordered_figures.extend(list(sorted(figs_order, key=figs_order.__getitem__)))

    return ordered_figures

def assign_data(styles, displayed_figs, data):
    """
       Add 'xdata', 'ydata' and 'legend_data' to 'lines_styles'
       and return 'figs_ids' of figures that do not contain
       and data.
     """

    # build a list of figs that actually contain some data
    nonempty_figs = []

    if not data.has_data('lines'):
        print('No data is provided. Nothing to display.')
        return nonempty_figs

    lines_ids = styles['lines_order']
    plots_keep = styles['plots_keep']
    plots_remove = styles['plots_remove']
    lines_styles = styles['lines']

    for line_id in lines_ids:

        line_data = data.get_linedata(line_id, not_found={})

        if not has_data(line_data): continue

        line_style = lines_styles[line_id]

        for fig_id in displayed_figs:
            # we skip other 'h' and 'u' data, as it would be mess
            if (fig_id in ['h', 'u']) and (not line_id == 'computed'):
                continue

            if not fig_id in line_data: continue

            line_value = line_data[fig_id]

            (xdata, ydata) = (line_value[0], line_value[1])

            if (not has_data(xdata)) or (not has_data(ydata)): continue

            # This fig surely will have some data
            if not fig_id in nonempty_figs:
                nonempty_figs.append(fig_id)

            if not fig_id in line_style:
                line_style[fig_id] = {}

            line_fig_style = line_style[fig_id]

            # Process special 'legend_data' (e.g. for 'h' and 'u' it's times)
            # as it needs to be converted and filtered
            if len(line_value) > 2:
                legend_data = line_value[2]
                if (type(legend_data) is np.ndarray):
                    legend_data = ['% 6.1f' % (ti/60.) for ti in legend_data]
            else:
                legend_data = None

            if ((not line_id in ('measured', 'original'))
                and ((fig_id in plots_keep) or (fig_id in plots_remove))):

                # filter computed data
                if fig_id in ('h', 'u'):
                    filter_size = xdata.shape[1]
                else:
                    filter_size = np.alen(xdata) # can be (n, ) ~ (1,n) ~ (n,1)

                if fig_id in plots_keep:
                    filter_idxs = np.zeros((filter_size, ), dtype=bool)
                    filter_indices(filter_idxs, plots_keep[fig_id], True)
                else:
                    filter_ones = np.zeros((filter_size, ), dtype=bool)

                if fig_id in plots_remove:
                    filter_indices(filter_idxs, plots_remove[fig_id], False)

                if fig_id in ('h', 'u'):
                    xdata = xdata[:, filter_idxs]
                    ydata = ydata[:, filter_idxs]
                else:
                    xdata = xdata[filter_idxs]
                    ydata = ydata[filter_idxs]

                if not legend_data is None:
                    if type(legend_data) is np.ndarray:
                        legend_data = legend_data[filter_idxs]

                    line_fig_style['legend_data'] = legend_data

            if fig_id in ('theta', 'relsat'):
                # swtich xdata <-> ydata
                tmp = xdata; xdata = ydata; ydata = tmp

            line_fig_style['xdata'] = xdata
            line_fig_style['ydata'] = ydata

    return nonempty_figs

def get_shown_figs_ids(figures_styles):
    # Keep only figures with 'show' == True
    figs_ids = [fig_id for fig_id in figures_styles.keys()
                  if get_figure_option(figures_styles, fig_id, 'show')]

    return figs_ids

class DPlots():
    """ Class for displaying data. User styles are applied. """

    def __init__(self, experiment_info):
        self._experiment_info = experiment_info

        # Default values (merged with plotstyles inifiles)
        display_options = {'separate_figures': False, 'show_figures': True,
                           'save_figures': True, 'matplotlib_backend': None}

        styles = {'options': DISPLAY_OPTIONS, 'lines': LINESTYLES_DEFAULT,
                  'figures': set_default_units(FIGURES_DEFAULTS),
                  'plots_keep': {}, 'plots_remove': {}, 'params_ref': {},
                  'lines_order': ()}

        # Read user plotysles inifiles
        plotstyles_filenames = get_filenames(experiment_info)

        user_styles = {}
        if plotstyles_filenames:
            for fname in plotstyles_filenames:
                deep_dictupdate(user_styles, read_plotstyles_file(fname))

        update_styles(styles, user_styles)

        self._styles    = styles

        display_options = styles['options']

        if 'matplotlib_backend' in display_options:
            matplotlib_backend = display_options['matplotlib_backend']
        if matplotlib_backend: # chosen other backend than the default one
            matplotlib.use(matplotlib_backend)

    def get_references(self):
        return self._styles['params_ref']

    def display(self, data, fignum=1):
        """ Display the figures and/or write them to files. """

        styles = self._styles

        display_options = styles['options']

        show_figures     = display_options['show_figures']
        save_figures     = display_options['save_figures']
        separate_figures = display_options['separate_figures']

        if not (save_figures or show_figures):
            return

        figs_styles = styles['figures']
        figs_ids = assign_data(styles, get_shown_figs_ids(figs_styles), data)

        if not figs_ids:
            return

        figs_ids = order_figures(figs_styles, figs_ids)

        lines_styles = styles['lines']
        lines_ids = styles['lines_order']

        print_status(data)

        import matplotlib.pyplot as plt

        if save_figures:
            experiment_info = self._experiment_info
            if experiment_info['mask']:
                figs_dir_type = 'mask'
            else:
                figs_dir_type = 'data'

            save_dir = get_directories('figs', figs_dir_type, experiment_info)

            if not path.exists(save_dir):
                makedirs(save_dir)

        if separate_figures:
            images_per_figure = 1
        else:
            images_per_figure = 6

        fignum -= 1
        img_num = 2^20 # high initialization, so that first fig is created

        for fig_id in figs_ids:
            # resolve figure and subplot
            if img_num > images_per_figure:
                img_num = 1
                fignum += 1

                if separate_figures:
                    plt.figure(fignum)
                else:
                    plt.figure(fignum, figsize=(16, 8.5))
                    plt.subplots_adjust(wspace=0.15, left=0.06, right=0.85)

            if not separate_figures:
                plt.subplot(3, 2, img_num)

            # plot the lines
            legend_labels = []

            xunit = get_figure_option(figs_styles, fig_id, 'xunit')
            yunit = get_figure_option(figs_styles, fig_id, 'yunit')

            for line_id in lines_ids:
                xdata = get_line_option(lines_styles, line_id, 'xdata', fig_id)
                ydata = get_line_option(lines_styles, line_id, 'ydata', fig_id)

                if (xdata is None) or (ydata is None):
                    continue

                plot_style = get_line_option(lines_styles, line_id, 'lineopt',
                                             fig_id)
                width = get_line_option(lines_styles, line_id, 'width', fig_id)
                symbolsize = get_line_option(lines_styles, line_id,
                                             'symbolsize', fig_id)
                xcoef = get_unit_coef(xunit)
                ycoef = get_unit_coef(yunit)
                if not xcoef == 1.0:
                    xdata = xcoef * np.asarray(xdata, dtype=float)
                if not ycoef == 1.0:
                    xdata = ycoef * np.asarray(ydata, dtype=float)
                ls = get_line_option(lines_styles, line_id, 'ls', fig_id)
                if ls and len(xdata.shape) > 1:
                    for ind in range(len(xdata[0])):
                        entryx = xdata[:, ind]
                        entryy = ydata[:, ind]
                        if symbolsize:
                            plt.plot(entryx, entryy, ls[ind%len(ls)],
                             linewidth=width, markersize=symbolsize)
                        else:
                            plt.plot(entryx, entryy, ls[ind%len(ls)],
                             linewidth=width)
                else:
                    if symbolsize:
                        plt.plot(xdata, ydata, plot_style, linewidth=width,
                                 markersize=symbolsize)
                    else:
                        plt.plot(xdata, ydata, plot_style, linewidth=width)

                # Extend the legend labels
                legend_label = get_line_option(lines_styles, line_id,
                                               'legend_data', fig_id)
                if legend_label is None:
                    legend_label =  get_line_option(lines_styles, line_id,
                                                    'label', fig_id, line_id)
                if np.isscalar(legend_label):
                    legend_labels.append(legend_label)
                else:
                    legend_labels.extend(legend_label)

            xlabel = get_figure_option(figs_styles, fig_id, 'xlabel')
            ylabel = get_figure_option(figs_styles, fig_id, 'ylabel')

            plt.xlabel(xlabel.format(xunit))
            plt.ylabel(ylabel.format(yunit))

            xscale = get_figure_option(figs_styles, fig_id, 'xscale')
            yscale = get_figure_option(figs_styles, fig_id, 'yscale')
            if xscale: plt.xscale(xscale)
            if yscale: plt.yscale(yscale)

            xmin = get_figure_option(figs_styles, fig_id, 'xmin')
            ymin = get_figure_option(figs_styles, fig_id, 'ymin')
            xmax = get_figure_option(figs_styles, fig_id, 'xmax')
            ymax = get_figure_option(figs_styles, fig_id, 'ymax')
            if np.isscalar(xmin): plt.xlim(xmin=xmin)
            if np.isscalar(xmax): plt.xlim(xmax=xmax)
            if np.isscalar(ymin): plt.ylim(ymin=ymin)
            if np.isscalar(ymax): plt.ylim(ymax=ymax)

            show_legend = get_figure_option(figs_styles, fig_id, 'show_legend')
            if show_legend:
                legend_loc   = get_figure_option(figs_styles, fig_id, 'legend_loc')
                legend_title = get_figure_option(figs_styles, fig_id, 'legend_title')
                legend_bbox  = get_figure_option(figs_styles, fig_id, 'legend_bbox')
                plt.legend(legend_labels, borderaxespad=0.0,
                           prop={'family': 'monospace'}, loc=legend_loc,
                           title=legend_title, bbox_to_anchor=legend_bbox)

            if save_figures and (img_num == images_per_figure):
                if separate_figures: img_suffix = fig_id
                else: img_suffix = str(fignum)

                plt.savefig(save_dir + 'image-' + img_suffix, dpi=300)

            img_num += 1

        if save_figures and (img_num < images_per_figure):
            plt.savefig(save_dir + 'image-' + str(fignum), dpi=300)

        if show_figures:
            try:
                plt.show(block=False)
            except: # Older matplotlib compatibility
                plt.ion()
                plt.show()

            if sys.version_info[0] == 2:
                #python 2.7
                raw_input('Press ENTER to continue...')
            else:
                input('Press ENTER to continue...')

            plt.close('all')

################################################################################
#                         Module's provided functions                          #
################################################################################

def show_results(experiment_info,
                 model=None, inv_params=None, cov=None):
    """
      Run the computation, store the results, load custom display styles
      and previously stored values and show/store the final figures.
    """

    data = DataStorage(experiment_info)

    if model is None:
        if not data.load():
            print('      (Was computation already run?)'
                  '\nINFO: Nothing to display. Exiting.')
            exit(0)

        save_data = False
    else:
        # Inverse problems set to compute data needed for running inverse;
        # Disable this to determine all data in the direct run
        model.measurements.run_inverse_problem_p(False)

        data.store_measurements(model.measurements, model)
        data.store_computation(model, model.measurements)

        if not inv_params is None: data.store('inv_params', inv_params)
        if not cov is None: data.store('cov', cov)

        savedir = get_directories('figs', 'mask', experiment_info)
        filename = savedir + '/' + 'results.txt'

        if not path.exists(savedir):
            makedirs(savedir)

        print_status(data, filename)

        save_data = True

    dplots = DPlots(experiment_info)

    if data.store_references(dplots.get_references(), model) or save_data:
        data.save()

    dplots.display(data)

################################################################################
#                              Unused functions                                #
################################################################################

def nd2strlist(nd):
    result = []
    for value in nd:
        result.append(str(value))
    return result

def display_table(t_measured=None, t_computed=None,
                  wl_out1_measured=None, wl_out1_computed=None,
                  gc1_measured=None, gc1_computed=None,
                  rm1_measured=None, rm1_computed=None,
                  l0_measured=None, l1_measured=None, l1_computed=None,
                  fignum=10):
    import matplotlib.pyplot as plt
    min_value = 1.0e-10
    assure = lambda v: max(v, min_value)
    format_row = (lambda format_str, data_row:
                  [format_str % assure(value) for value in data_row])

    disp_t      = (not t_measured is None) and (not t_computed is None)
    disp_wl_out = ((not wl_out1_measured is None)
                   and (not wl_out1_computed is None))
    disp_gc     = (not gc1_measured is None) and (not gc1_computed is None)
    disp_rm     = (not rm1_measured is None) and (not rm1_computed is None)
    disp_l      = (not l1_measured is None) and (not l1_computed is None)

    disp_p = (disp_t, disp_wl_out, disp_gc, disp_rm, disp_l)
    disp_labels = ('t', 'wl_out', 'gc', 'rm', 'l')

    subplots = sum([int(value) for value in disp_p])
    print('sb', subplots)

    colLabels = ['#%i' % (i+1) for i in range(len(t_measured))]
    print(wl_out1_measured, wl_out1_computed, colLabels)

    plt.figure(fignum, figsize=(16, 8.5))

    subplt_num = 1

    for (disp, label) in zip (disp_p, disp_labels):
        if not disp: continue

        data = []

        if label == 't':
            rowLabels = ['Duration measured', 'Duration computed', 'Error (%)']
            data_measured = t_measured
            data_computed = t_computed
        elif label == 'wl_out':
            rowLabels = ['Outflow measured', 'Outflow computed', 'Error (%)']
            data_measured = wl_out1_measured
            data_computed = wl_out1_computed
        elif label == 'gc':
            rowLabels = ['GC measured', 'GC computed', 'Error (%)']
            data_measured = gc1_measured
            data_computed = gc1_computed
        elif label == 'rm':
            rowLabels = ['RM measured', 'RM computed', 'Error (%)']
            data_measured = rm1_measured
            data_computed = rm1_computed
        elif label == 'l':
            rowLabels = ['L1 measured', 'L1 computed', 'Error (%)']
            if not l0_measured is None:
                rowLabels.insert(0, 'L0 initial')
            data_measured = l1_measured
            data_computed = l1_computed

        plt.subplot(subplots, 1, subplt_num)
        #plt.axes(frameon=False)
        plt.axis('off')

        data.append(format_row('% 9.6f', data_measured))
        data.append(format_row('% 9.6f', data_computed))
        data.append(format_row('% 5.2f',
                               [((wo_m - wo_c) / wo_m * 100) for (wo_m, wo_c)
                                in zip(data_measured, data_computed)]))

        subplt_num = subplt_num + 1
        print(data)

        plt.table(cellText=data,
                  rowLabels=rowLabels,
                  colLabels=colLabels,
                  loc='top')
    try:
        plt.show(block=False)
    except:
        plt.ion()
        plt.show()
    input('Press ENTER to continue...')
