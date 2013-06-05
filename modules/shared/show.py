from __future__ import print_function, division

import numpy as np, matplotlib, pickle
import sys
from const import PLOTSTYLE_ININAME, DUMP_DATA_VERSION, DUMP_DATA_FILENAME
from os import makedirs, path
from shared import get_directories, parse_value
from config import ModulesManager, load_model
from collections import OrderedDict
from modules.shared.functions import has_data, compare_data
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
              'none': ('', )}
DEFAULT_UNITS = {'length': 'cm', 'time': 'min', 'pressure': 'Pa',
                 'force_kgp': 'gf', 'weight': 'g', 'velocity': 'cm/s',
                 'none': ''}

dg_label_time = "Time [{}]"
dg_label_length = "Sample length $L$ [{}]"
dg_unit_time_length = ('time', 'length')
DG_AXES_LABELS = \
  OrderedDict({'h': ((dg_label_length, "Piezometric head $h$ [{}]"),
                     ('length', 'length')),
               'u': ((dg_label_length, "Effective saturation $S_e${}"),
                     ('length', 'none')),
               'MO': ((dg_label_time, "Expelled water [{}]"),
                      dg_unit_time_length),
               'MI': ((dg_label_time, "Inflow water [{}]"),
                      dg_unit_time_length),
               'GC': ((dg_label_time, "Gravitational center [{}]"),
                      dg_unit_time_length),
               'RM': ((dg_label_time, "Rotational momentum [kg.m.s$^{-1}$]{}"),
                      ('time', 'none')),
               's1': ((dg_label_time, "Interface s1 [{}]"),
                      dg_unit_time_length),
               's2': ((dg_label_time, "Interface s2 [{}]"),
                      dg_unit_time_length),
               'WM': ((dg_label_time, "Water mass balance [{}]"),
                      dg_unit_time_length),
               'WM_in_tube': ((dg_label_time, "Water mass in tube [{}]"),
                              dg_unit_time_length),
               'theta': (("Water content $\\theta${}", "Pressure $p$ [{}]"),
                         ('none', 'pressure')),
               'relsat': (("Effective saturation $S_e${}", "Hydraulic head $h$ [{}]"),
                          ('none', 'length')),
               'K':  (("Water content $\\theta${}",
                       "Hydraulic conductivity $K(\\theta)$ [{}]"),
                      ('none', 'velocity')),
               'gF_MO': ((dg_label_time, "Force of expelled water [{}]"),
                         ('time', 'force_kgp')),
               'gF_MT': ((dg_label_time, "Force of water in tube [{}]"),
                         ('time', 'force_kgp')),
               'dgF_MO': ((dg_label_time, "Force difference of expelled water [{}]"),
                          ('time', 'force_kgp')),
               'dgF_MT': ((dg_label_time, "Force difference of water in tube [{}]"),
                          ('time', 'force_kgp'))})

DG_PAIRS = (('h', 'u'), ('MI', 'MO'), ('GC', 'RM'), ('gF_MT', 'gF_MO'),
            ('dgF_MT', 'dgF_MO'), ('s1', 's2'), ('theta', 'relsat'))

FIGURES_IDS = list(DG_AXES_LABELS.keys())

def get_unit_coef(unit_base):
    unit = unit_base.lower()
    # units used for computation are: cm, s, pa, gf and "no units"
    if unit in ['cm', 's', 'pa', 'g', 'gf', 'cm/s', '']: coef = 1.0
    elif unit == 'mm': coef = 10.
    elif unit == 'min': coef = 1./60.
    elif unit == 'h': coef = 1./3600.
    elif unit in ('m', 'm/s'): coef = 0.01
    elif unit in ['kpa', 'kg', 'kgp']: coef = 0.001
    else:
        print('Unknown unit:', unit_base, '\nKnown units are only:', DATA_UNITS)
        exit(1)
    return coef

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
        return ((data_type in self._data)
                and (not self._data[data_type] is None))

    def store(self, name, value):
        self._data[name] = value

    def get(self, key, not_found=None):
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

                if hasattr(model, 'p') or hasattr(model, 'pc'):
                    if hasattr(model, 'p'):
                        p_meas = np.asarray(model.p)
                    else:
                        p_meas = np.asarray(model.pc)

                    (p_user, theta_user) = \
                      SC.retention_curve(theta_s, model.density, model.g,
                                         theta_r=theta_r, p=p_meas)

                    data['theta'] = (p, theta, theta_user)
                else:
                    data['theta'] = (p, theta)

                K = SC.conductivity_curve(model.ks, theta_s,
                                          theta_r=theta_r, g=model.g,
                                          rho=model.density)
                data['K'] = K

            self._data['lines'][ID] = data
            self.store('experiment_info', model.experiment_info)

        return flag

    def store_references(self, user_references, model=None):
        stored_references = self.get('references', not_found={})

        if type(user_references) in (list, tuple):
            print("DEPRECATION ERROR: format of 'params_ref' has been changed."
                  "Please update to new format.")
            return False

        if user_references == stored_references:
            return False   # nothing needs to be re-computed

        if model is None:
            model = load_model(self.get('experiment_info'), validate=True)

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
            model.set_parameters(ref_params)

            flag = self.store_computation(model, model.measurements, ID=ref_id)

            if not flag:
                print('Reference parameters: ', ref_params)

            model.set_parameters(backup_params) # restore

        self.store('references', stored_references)

        return True

    def store_measurements(self, measurements, model=None):
        for mtype in ('original', 'measured'):
            m = {}

            untransformed = (mtype == 'original')
            for (name, xvalue, yvalue) \
                in measurements.iterate_meas_measurements(untransformed, model):

                m[name] = (xvalue, yvalue)

                self._data['lines'][mtype] = m

    def get_linedata(self, line_id, not_found=None):
        data = self._data['lines']

        if line_id in data:
            return data[line_id]
        else:
            return not_found

    def save(self):
        if not self._data:
            print('No data was stored. Nothing to be saved. Skipping saving...')
            return

        savedir = get_directories('figs', 'mask', self.get('experiment_info'))
        if not path.exists(savedir):
            makedirs(savedir)

        with open(savedir + DUMP_DATA_FILENAME, 'wb') as fout:
            pickle.dump(self._data, fout, DUMP_DATA_VERSION)

    def load(self):
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

def mk_status_item(data_id, data_computed, data_measured = []):
    return {'id': data_id, 'data': (data_computed, data_measured)}

def display_status(data_plots=None, stream=None):

    if not data_plots: return

    for plot in data_plots:
        plot_id = plot['id']

        data_items_nr = len(plot['data'])
        if (data_items_nr == 0 ) or (not has_data(plot['data'][0])):
            continue

        value_computed = np.asarray(plot['data'][0])

        if (data_items_nr == 1 ) or (not has_data(plot['data'][1])):
            value_measured =  None
        else:
            value_measured = np.asarray(plot['data'][1])

        compare_data(plot_id, value_computed, value_measured, stream)

def print_status(data, filename=None):
    computed     = data.get_linedata('computed')
    measurements = data.get_linedata('measured')

    if not computed: return

    if filename is None:
        stream = None
    else:
        stream = open(filename, 'w')

    status_items = []

    if not measurements: return

    try:
        # compare measured vs. computed data
        for (key, m_data) in measurements.items():
            if m_data[1] is None: continue

            if key in computed:
                if key == 'theta':
                    c_value = computed[key][2]
                else:
                    c_value = computed[key][1][1:]
                m_value = m_data[1]
            else:
                continue

            if c_value is not None:
                status_items.append(mk_status_item(key, c_value, m_value))

        if status_items:
            display_status(status_items, stream)

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
            print('WARNING: Errors during parsing occured, only default values '
                  'will be used.')
            result = {} # whatever was read so far is discared
            break

    return result

def mk_figurestyles(fig_id):
    """ Add default styles for figure with given 'fig_id' """

    figure_styles = dict.fromkeys(('xlabel', 'ylabel',
                                   'xscale', 'yscale', 'xunit', 'yunit',
                                   'xmin', 'ymin', 'xmax', 'ymax',
                                   'show', 'ls', 'show_legend',
                                   'legend_title', 'legend_bbox', 'legend_loc'))

    # Default values
    if fig_id in DG_AXES_LABELS:
        figure_styles['xlabel'] = DG_AXES_LABELS[fig_id][0][0]
        figure_styles['ylabel'] = DG_AXES_LABELS[fig_id][0][1]

        unit_type = DG_AXES_LABELS[fig_id][1][0]
        figure_styles['xunit'] = DEFAULT_UNITS[unit_type]

        unit_type = DG_AXES_LABELS[fig_id][1][1]
        figure_styles['yunit'] = DEFAULT_UNITS[unit_type]

        figure_styles['legend_loc'] = 4


    if fig_id in ['h', 'u']:
        figure_styles['legend_title'] = 'Time [min]'

    elif fig_id == 'theta':
        figure_styles['yscale'] = 'log'
        figure_styles['legend_loc'] = 1

    elif fig_id == 'K':
        figure_styles['yscale'] = 'log'

    return figure_styles

def figures_styles_post_update(figures_styles, display_options):
    """
      This function is called after merging the default plotting styles
      with styles from plotstyles inifiles. Styles can be adapted based
      on the information in the inifiles.
    """

    h_styles = figures_styles['h']

    if ((h_styles['show_legend'] is None)
        and (not display_options['separate_figures'])):
        h_styles['show_legend'] = False
    elif (h_styles['legend_loc'] is None) and (h_styles['legend_bbox']):
        h_styles['legend_bbox'] = (1.01, 1.)
        h_styles['legend_loc'] = 2

    u_styles = figures_styles['u']
    if (u_styles['legend_loc'] is None) and (u_styles['legend_bbox']):
        u_styles['legend_bbox'] = (1.01, 1.)
        u_styles['legend_loc'] = 2

def linestyles_post_update(styles):
    """
      Lines are specified in 'params_ref' variable of the plotstyles
      inifiles, hence as we don't know them apriori (of course except for
      lines 'measured and 'computed') and hence only post-update method is
      available - all the default values for all lines are therefore specified
      here.
      The values are returned in the final (correct) format (which may be
      different from the values read from the inifiles).
      Correct lines order based on the 'order' option for the line is
      determined here too.
    """

    lines_ids = (['original', 'measured', 'computed']
                 + list(styles['params_ref'].keys()))
    lines_styles = styles['lines']
    lines_order  = {}

    for line_id in lines_ids:
        # Determine user supplied line styles
        if not line_id in lines_styles:
            lines_styles[line_id] = {}

        line_styles = lines_styles[line_id]

        # Process line label
        if 'label' in line_styles:
            label = line_styles['label']
            del line_styles['label']
        else:
            label = line_id

        # Process line width
        if 'width' in line_styles:
            width = line_styles['width']
        else:
            width = 1

        if 'symbolsize' in line_styles:
            symbolsize = line_styles['symbolsize']
        else:
            symbolsize = None
        # Process line order
        if 'order' in line_styles:
            lines_order[line_id] = line_styles['order']
            del line_styles['order']
        else:
            lines_order[line_id] = 999

        # Default lineopt
        if line_id == 'measured':
            default_lineopt = 'x'
        else:
            default_lineopt = '.'

        # Set correct format for line_styles and supply (default) missing values
        for fig_id in FIGURES_IDS:
            if fig_id in line_styles:
                lineopt = line_styles[fig_id]
            elif fig_id in ('h', 'u'):
                lineopt = '-'
            elif (fig_id == 'theta') and (not line_id == 'measured'):
                lineopt = '-'
            else:
                lineopt = default_lineopt

            line_styles[fig_id] = {'lineopt': lineopt, 'label': label, 'width': width,
                            'symbolsize': symbolsize}

    ordered_lines = list(sorted(lines_order, key=lines_order.__getitem__))
    styles['lines_order'] = ordered_lines

def order_figures(figures):
    """
      Sort figures according to prefered order. Some figures are prefered
      to be next to each other (see DG_PAIRS variable).
    """

    # Remove figures with no data and not displayed
    for fig_id in list(figures.keys()):
        if ((not has_data(figures[fig_id]['data']))
            or (figures[fig_id]['styles']['show'] == False)):
            del figures[fig_id]

    # Order remaining figures
    ordered_figures = []

    # first order pairs (and queue singles from pairs)
    for (i1_id, i2_id) in DG_PAIRS:
        if (i1_id in figures) and (i2_id in figures):
            ordered_figures.append(figures[i1_id])
            ordered_figures.append(figures[i2_id])

            del figures[i1_id]
            del figures[i2_id]

    for fig_id in FIGURES_IDS:
        if fig_id in figures:
            ordered_figures.append(figures[fig_id])

    return ordered_figures

def mk_figures(data, styles):
    """ Create figures based on the supplied data and styles information. """

    if not data.has_data('lines'):
        print('No data is provided. Nothing to display.')
        return None

    figures = {fig_id: {'id': fig_id, 'data': [],
                        'styles': styles['figures'][fig_id]}
                for fig_id in FIGURES_IDS}

    lines_ids = styles['lines_order']

    for line_id in lines_ids:
        line_data = data.get_linedata(line_id, not_found={})
        lines_styles = styles['lines'][line_id]

        for (fig_id, line_value) in line_data.items():
            # Filter out supplied data not recognized as valid
            if (not fig_id in FIGURES_IDS) or (not has_data(line_data)):
                continue

            # we skip other 'h' and 'u' data, as it would be mess
            if (fig_id in ['h', 'u']) and (not line_id == 'computed'):
                continue

            (xdata, ydata) = (line_value[0], line_value[1])

            if not (has_data(xdata) and has_data(ydata)): continue

            line_style = lines_styles[fig_id]

            if ((fig_id in ['h', 'u']) and (len(line_value) > 2)
                and has_data(line_value[2])):
                line_style['label'] = \
                  ['% 6d' % (ti/60.) for ti in line_value[2]]

            if fig_id == 'theta':
                item = {'xdata': ydata, 'ydata': xdata}
            else:
                item = {'xdata': xdata, 'ydata': ydata}

            item.update(line_style)

            figures[fig_id]['data'].append(item)

    return order_figures(figures)

class DPlots():
    """ Class for displaying data. User styles are applied. """

    def __init__(self, experiment_info):
        self._experiment_info = experiment_info

        # User styles (read from plotstyles inifiles)
        display_options = {'separate_figures': False, 'show_figures': True,
                           'save_figures': True, 'matplotlib_backend': None}
        figures_styles = {fig_id: mk_figurestyles(fig_id)
                          for fig_id in FIGURES_IDS}

        styles = {'options': display_options, 'figures': figures_styles,
                  'lines': {}, 'params_ref': {}, 'lines_order': ()}

        plotstyles_filenames = get_filenames(experiment_info)

        if plotstyles_filenames:
            for fname in plotstyles_filenames:
                deep_dictupdate(styles, read_plotstyles_file(fname))

        figures_styles_post_update(figures_styles, display_options)
        linestyles_post_update(styles)

        self._styles    = styles

        # Global displaying options
        self._display_options = display_options

        if 'matplotlib_backend' in display_options:
            matplotlib_backend = display_options['matplotlib_backend']
        if matplotlib_backend: # chosen other backend than the default one
            matplotlib.use(matplotlib_backend)

    def get_references(self):
        return self._styles['params_ref']

    def display(self, data, fignum=1):
        """ Display the figures and/or write them to files. """

        display_options = self._display_options

        show_figures     = display_options['show_figures']
        save_figures     = display_options['save_figures']
        separate_figures = display_options['separate_figures']

        if not (save_figures or show_figures):
            return

        dplots = mk_figures(data, self._styles)

        if not dplots:
            return

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

        for figure in dplots:
            if not figure['data']: continue
            fig_id = figure['id']

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

            # plot the supplied data
            figure_styles = figure['styles']
            (xunit, yunit) = (figure_styles['xunit'], figure_styles['yunit'])

            plot_labels = []
            for figure_data in figure['data']:
                xdata = figure_data['xdata']
                ydata = figure_data['ydata']
                line_label = figure_data['label']
                plot_style = figure_data['lineopt']
                width = figure_data['width']
                symbolsize = figure_data['symbolsize']

                xcoef = get_unit_coef(xunit)
                ycoef = get_unit_coef(yunit)
                if not xcoef == 1.0:
                    xdata = xcoef * np.asarray(xdata, dtype=float)
                if not ycoef == 1.0:
                    xdata = ycoef * np.asarray(ydata, dtype=float)
                if figure_styles['ls'] and len(xdata.shape) > 1:
                    for ind in range(len(xdata[0])):
                        entryx = xdata[:, ind]
                        entryy = ydata[:, ind]
                        if symbolsize:
                            plt.plot(entryx, entryy, 
                             figure_styles['ls'][ind%len(figure_styles['ls'])],
                             linewidth=width, markersize=symbolsize)
                        else:
                            plt.plot(entryx, entryy, 
                             figure_styles['ls'][ind%len(figure_styles['ls'])],
                             linewidth=width)
                else:
                    if symbolsize:
                        plt.plot(xdata, ydata, plot_style, linewidth=width,
                                    markersize=symbolsize)
                    else:
                        plt.plot(xdata, ydata, plot_style, linewidth=width)
                if type(line_label) == str:
                    plot_labels.append(line_label)
                else:
                    plot_labels.extend(line_label)

            xlabel = figure_styles['xlabel']
            ylabel = figure_styles['ylabel']

            plt.xlabel(xlabel.format(xunit))
            plt.ylabel(ylabel.format(yunit))
            if figure_styles['xscale']:
                plt.xscale(figure_styles['xscale'])
            if figure_styles['yscale']:
                plt.yscale(figure_styles['yscale'])

            if figure_styles['xmin']:
                plt.xlim(xmin=figure_styles['xmin'])
            if figure_styles['xmax']:
                plt.xlim(xmax=figure_styles['xmax'])
            if figure_styles['ymin']:
                plt.ylim(ymin=figure_styles['ymin'])
            if figure_styles['ymax']:
                plt.ylim(ymax=figure_styles['ymax'])

            show_legend = figure_styles['show_legend']
            if show_legend is None:
                if (len(figure['data']) > 1) or (np.ndim(ydata) > 1):
                    show_legend = True
            if show_legend:
                plt.legend(plot_labels, borderaxespad=0.0,
                           prop={'family': 'monospace'},
                           loc=figure_styles['legend_loc'],
                           title=figure_styles['legend_title'],
                           bbox_to_anchor=figure_styles['legend_bbox'])

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
