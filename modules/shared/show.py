from __future__ import print_function, division

import matplotlib
import sys
import numpy as np
from const import FIGS_DIR, PLOTSTYLE_ININAME, DUMP_DATA_VERSION
from os import makedirs, path
from shared import get_directories, parse_value
from config import ModulesManager, load_model, DataStorage
from modules.shared.functions import has_data, compare_data
try:
    import ConfigParser as configparser
except:
    import configparser

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
DG_AXES_LABELS = {'h': ((dg_label_length, "Piezometric head $h$ [{}]"),
                        ('length', 'length')),
                  'u': ((dg_label_length, "Relative saturation $u${}"),
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
                  'relsat': (("Relative saturation $u${}", "Hydraulic head $h$ [{}]"),
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
                           ('time', 'force_kgp'))}
DG_PAIRS = (('h', 'u'), ('MI', 'MO'), ('GC', 'RM'), ('gF_MT', 'gF_MO'),
            ('dgF_MT', 'dgF_MO'), ('s1', 's2'), ('theta', 'relsat'))

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
    if filename is None:
        stream = None
    else:
        stream = open(filename, 'w')

    status_items = []

    measurements = data.get_linedata('measured')
    computed     = data.get_linedata('computed')

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
        cov = data.get_value('cov')
        if not cov is None:
            print('\nCov:\n', cov, file=stream)
            print('Deviation:\n', np.sqrt(np.diag(cov)), file=stream)

        params = data.get_value('inv_params')
        if not params is None:
            print('\nOptimal parameters found:', file=stream)
            for (name, value) in params.items():
                if name == 'ks':
                    print('  Ks [cm/s]: {: .8g}'.format(value), file=stream)
                else:
                    print('  {:9}: {: .8g}'.format(name, value), file=stream)

    finally:
        if not filename is None: stream.close()

# ResultData: hold the data of the computation
# Structure: {'lines': lines_structure, 'inv_params': inv_params, 'cov': cov}
# where:
# lines_structure: dictionary of types: {line_id: line_data} where:
#     line_id the ID of the line
#     line_data a dict of types {data_type: (xdata, ydata)}, where:
#          data_type is (by default) one of ['h', 'u', 'GC', 'RM', 'WM', 's1','s2'],
#          xdata (ydata)is the x-axis (y-axis) coordinate
class ResultsData():
    def __init__(self):
        self._data = {'lines': {}}
        self._modman = None

    def has_data(self, data_type):
        return ((data_type in self._data)
                and (not self._data[data_type] is None))

    def store_value(self, name, value):
        self._data[name] = value

    def get_value(self, name, not_found=None):
        if name in self._data:
            return self._data[name]
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
                from modules.shared.saturation_curve import (retention_curve,
                                                             conductivity_curve)

                if hasattr(model, 'theta_s'): theta_s = model.theta_s
                else: theta_s = model.porosity

                if hasattr(model, 'theta_r'): theta_r = model.theta_r
                else: theta_r = 0.0

                (p, theta) = retention_curve(model.SC,
                                             theta_s, model.density, model.g,
                                             theta_r=theta_r)
                if hasattr(model, 'p'):
                    (p_user, theta_user) = \
                      retention_curve(model.SC, theta_s,
                                      model.density, model.g, theta_r=theta_r,
                                      p=None, h=model.h)
                    data['theta'] = (p, theta, theta_user)
                else:
                    data['theta'] = (p, theta)

                data['K'] = conductivity_curve(model.SC, model.ks, theta_s,
                                               theta_r=theta_r, g=model.g,
                                               rho=model.density)

            self._data['lines'][ID] = data
            self.store_value('experiment_info', model.experiment_info)

        return flag

    def store_references(self, user_references, model=None):
        stored_references = self.get_value('references', not_found={})
        data = self._data['lines']

        if type(user_references) in (list, tuple):
            print("DEPRECATION ERROR: format of 'params_ref' has been changed."
                  "Please update to new format.")
            return False

        if user_references == stored_references:
            return False   # nothing needs to be re-computed

        if model is None:
            from config import load_model
            model = load_model(self.get_value('experiment_info'), validate=True)

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
                print('Reference parameters: ', ref)

            model.set_parameters(backup_params) # restore

        self.store_value('references', stored_references)

        return True

    def store_measurements(self, measurements):
        m = {}

        for (name, xvalue, yvalue) in measurements.iterate_meas_measurements():
            m[name] = (xvalue, yvalue)

        self._data['lines']['measured'] = m

    def get_linedatatypes(self):
        return DG_AXES_LABELS.keys()

    def get_linedata(self, line_id, not_found=None):
        data = self._data['lines']
        if line_id in data:
            return data[line_id]
        else:
            return not_found

    def iterate_lines(self):
        lines = self.get_value('lines')

        if lines is None:
            yield None
        else:
            for (line_id, line_data) in lines.items():
                yield (line_id, line_data)

    def dump(self):
        return self._data

    def load(self, value):
        if not value is None:
            self._data = value

################################################################################
#                         Process plotstyles.ini                               #
################################################################################
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
    """ Recursively update dictionary d1 with dictionary d2. """

    for (vkey, vvalue) in d2.items():
        if ((vkey in d1) and (type(d1[vkey]) == dict)
            and (type(vvalue) == dict)):
            deep_dictupdate(d1[vkey], vvalue)
        else:
            d1[vkey] = vvalue

    return d1

def read_plotstyles_file(filename):
    result = {}

    # read config files
    parser   = configparser.ConfigParser()
    try:
        read_files = parser.read(filename)
    except configparser.ParsingError as E:
        print('Plotstyles file ', filename,
              '\n could not be read because of the following error: ', E,
              '\nDefault values are used.')

        return result

    # Write data from parser to configuration
    for psection in parser.sections():
        section = psection.lower() # we ignore sections

    try:
        for option in parser.options(psection):
            raw_value = parser.get(psection, option).strip()

            result[option] = parse_value(raw_value, raise_error=True)
    except:
        print('WARNING: Errors during parsing occured, only default values '
              'will be used.')
        results = {}

    return result

def get_line_order(line_options, not_found=999):
    if 'order' in line_options:
        return line_options['order']
    else:
        return not_found

class PlotStyles():
    def __init__(self, experiment_info):
        # Load user styles (from plotstyles inifiles)
        plotstyles_filenames = get_filenames(experiment_info)

        if plotstyles_filenames:
            userstyles = read_plotstyles_file(plotstyles_filenames[0])

            for fname in plotstyles_filenames[1:]:
                deep_dictupdate(userstyles, read_plotstyles_file(fname))

        else:
            userstyles = {}

        self._userstyles = userstyles

        # Process global displaying options
        options = {'separate_figures': False, 'show_figures': True,
                   'save_figures': True, 'matplotlib_backend': None}
        user_opts = self.get_userstyles('options')
        if user_opts: options.update(user_opts)

        self._display_options = options
        self._figuresstyles = {}

    def get_userstyles(self, key):
        if key in self._userstyles:
            return self._userstyles[key]
        else:
            return None

    def get_display_options(self):
        return self._display_options

    def _mk_figuresstyles(self, dplot_id):
        dplot_styles = {tag: None for tag in ['xlabel', 'ylabel',
                                              'xscale', 'yscale',
                                              'xunit', 'yunit',
                                              'show_legend', 'legend_title',
                                              'legend_bbox', 'legend_loc',
                                              'show']}

        # Default values
        if dplot_id in DG_AXES_LABELS:
            dplot_styles['xlabel'] = DG_AXES_LABELS[dplot_id][0][0]
            dplot_styles['ylabel'] = DG_AXES_LABELS[dplot_id][0][1]

            unit_type = DG_AXES_LABELS[dplot_id][1][0]
            dplot_styles['xunit'] = DEFAULT_UNITS[unit_type]

            unit_type = DG_AXES_LABELS[dplot_id][1][1]
            dplot_styles['yunit'] = DEFAULT_UNITS[unit_type]

            dplot_styles['legend_loc'] = 4


        if dplot_id in ['h', 'u']:
            dplot_styles['legend_bbox'] = (1.01, 1.)
            dplot_styles['legend_loc'] = 2
            dplot_styles['legend_title'] = 'Time [min]'

        elif dplot_id == 'theta':
            dplot_styles['yscale'] = 'log'
            dplot_styles['legend_loc'] = 1

        elif dplot_id == 'K':
            dplot_styles['yscale'] = 'log'

        # User values specified in plotstyles file
        user_styles = self.get_userstyles('datasets')
        if user_styles and (dplot_id in user_styles):
            dplot_styles.update(user_styles[dplot_id])

        # Post-update
        if (dplot_id == 'h') and (dplot_styles['show_legend'] is None)
            and (not self.get_display_options()['separate_figures'])):
            dplot_styles['show_legend'] = False

        return dplot_styles

    def get_figuresstyles(self, dtype):
        if not dtype in self._figuresstyles:
            self._figuresstyles[dtype] = self._mk_figuresstyles(dtype)

        return self._figuresstyles[dtype]

    def get_linestyle(self, line_id, data_types):
        # set default values
        if line_id == 'measured':
            lineopt = 'x'
        else:
            lineopt = '.'

        line_styles = {dtype: lineopt for dtype in data_types}

        if 'h' in data_types: line_styles['h'] = '-'
        if 'u' in data_types: line_styles['u'] = '-'
        if ('theta' in data_types) and (not line_id == 'measured'):
             line_styles['theta'] = '-'

        # merge with user-specified values
        user_styles = self.get_userstyles('lines')
        if user_styles and (line_id in user_styles):
            line_styles.update(user_styles[line_id])

        return line_styles

################################################################################
#                             Data displaying                                  #
################################################################################
class DPlots():
    def __init__(self, experiment_info):
        self._experiment_info = experiment_info

        self._plotstyles = PlotStyles(experiment_info)

        self.fignum = 1

        display_options = self._plotstyles.get_display_options()
        if 'matplotlib_backend' in display_options:
            matplotlib_backend = display_options['matplotlib_backend']
        if matplotlib_backend: # chosen other backend than the default one
            matplotlib.use(matplotlib_backend)

    def get_references(self):
        return self._plotstyles.get_userstyles('params_ref')

    def get_lines_order(self):
        lines = self._plotstyles.get_userstyles('lines')
        refs  = self._plotstyles.get_userstyles('params_ref')

        valid_line_ids = ['measured', 'computed'] + list(refs.keys())

        # determine all IDs specified in 'lines' option and sort them by
        # their 'order' value
        lines_orders = [(line_id, get_line_order(line_data))
                        for (line_id, line_data) in lines.items()
                        if line_id in valid_line_ids]
        ordered_lines = sorted(lines_orders, key=lambda x: x[1])
        # ...and keep IDs only
        ordered_line_ids = [line_order[0] for line_order in ordered_lines]
        # also include IDs not specified in 'lines' options
        missing_line_ids = [line_id for line_id in valid_line_ids
                            if not line_id in ordered_line_ids]

        return ordered_line_ids + missing_line_ids

    def _mk_dplots(self, data):
        def _mk_dplots_bucket(data_types, plot_styles):
            dplots_bucket = \
              {dtype: {'id': dtype, 'data': [],
                       'styles': plot_styles.get_figuresstyles(dtype)}
                for dtype in data_types}

            return dplots_bucket

        def _add_plotline(line_id, line_data, data_types, plot_styles,
                          dplots_bucket):

            line_styles = plot_styles.get_linestyle(line_id, data_types)
            if 'label' in line_styles:
                label = line_styles['label']
            else:
                label = line_id

            for (data_type, data_value) in line_data.items():

                # Filter out supplied data not recognized as valid
                if not data_type in data_types: continue

                if not has_data(data_value): continue

                # we skip other 'h' and 'u' data, as it would be mess
                if (data_type in ['h', 'u']) and (not line_id == 'computed'):
                    continue

                (xdata, ydata) = (data_value[0], data_value[1])

                if not (has_data(xdata) and has_data(ydata)): continue

                if ((data_type in ['h', 'u']) and (len(data_value) > 2)
                    and has_data(data_value[2])):
                    ilabel = ['% 6d' % (ti/60.) for ti in data_value[2]]
                else:
                    ilabel = label

                if data_type == 'theta':
                    item = (ydata, xdata, ilabel, line_styles[data_type])
                else:
                    item = (xdata, ydata, ilabel, line_styles[data_type])
                dplots_bucket[data_type]['data'].append(item)

        def _filter_dplots(dplots_bucket):
            for name in list(dplots_bucket.keys()):
                if ((not has_data(dplots_bucket[name]['data']))
                    or (dplots_bucket[name]['styles']['show'] == False)):
                    del dplots_bucket[name]

            return dplots_bucket

        def _order_dplots(dplots_bucket):
            ordered_dplots = []
            single_dplots  = []

            # first order pairs (and queue singles from pairs)
            for (i1_id, i2_id) in DG_PAIRS:
                i1p = (i1_id in dplots_bucket)
                i2p = (i2_id in dplots_bucket)

                if i1p and i2p:
                    ordered_dplots.append(dplots_bucket[i1_id])
                    ordered_dplots.append(dplots_bucket[i2_id])
                elif i1p:
                    single_dplots.append(dplots_bucket[i1_id])
                elif i2p:
                    single_dplots.append(dplots_bucket[i2_id])

                if i1p: del dplots_bucket[i1_id]
                if i2p: del dplots_bucket[i2_id]

            # append singles
            ordered_dplots.extend(single_dplots)
            # append remaning singles that are not part of any pair
            for remaining_dplots in dplots_bucket.values():
                ordered_dplots.append(remaining_dplots)

            return ordered_dplots

        # function body
        if not data.has_data('lines'):
            print('No data is provided. Nothing to display.')
            return None

        plot_styles   = self._plotstyles
        data_types    = data.get_linedatatypes()
        dplots_bucket = _mk_dplots_bucket(data_types, plot_styles)

        ordered_line_ids = self.get_lines_order()

        for line_id in ordered_line_ids:
            _add_plotline(line_id, data.get_linedata(line_id), data_types,
                          plot_styles, dplots_bucket)

        ordered_dplots = _order_dplots(_filter_dplots(dplots_bucket))

        return ordered_dplots

    def display(self, data, fignum=1):

        display_options = self._plotstyles.get_display_options()

        show_figures     = display_options['show_figures']
        save_figures     = display_options['save_figures']
        separate_figures = display_options['separate_figures']

        if not (save_figures or show_figures):
            return

        dplots = self._mk_dplots(data)

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

        for dplot in dplots:
            if not dplot['data']: continue
            dplot_id = dplot['id']

            # resolve figure and subplot
            if img_num > images_per_figure:
                img_num = 1
                self.fignum += 1

                if separate_figures:
                    plt.figure(self.fignum)
                else:
                    plt.figure(self.fignum, figsize=(16, 8.5))
                    plt.subplots_adjust(wspace=0.15, left=0.06, right=0.85)

            if not separate_figures:
                plt.subplot(3,2,img_num)

            # plot the supplied data
            dplot_styles = dplot['styles']
            (xunit, yunit) = (dplot_styles['xunit'], dplot_styles['yunit'])

            plot_labels = []
            for dplot_data in dplot['data']:
                (xdata, ydata, line_label, plot_style) = dplot_data
                xcoef = get_unit_coef(xunit)
                ycoef = get_unit_coef(yunit)
                if not xcoef == 1.0:
                    xdata = xcoef * np.asarray(xdata, dtype=float)
                if not ycoef == 1.0:
                    xdata = ycoef * np.asarray(ydata, dtype=float)
                plt.plot(xdata, ydata, plot_style)
                if type(line_label) == str:
                    plot_labels.append(line_label)
                else:
                    plot_labels.extend(line_label)

            (xlabel, ylabel) = (dplot_styles['xlabel'], dplot_styles['ylabel'])
            plt.xlabel(xlabel.format(xunit))
            plt.ylabel(ylabel.format(yunit))
            if dplot_styles['xscale']:
                plt.xscale(dplot_styles['xscale'])
            if dplot_styles['yscale']:
                plt.yscale(dplot_styles['yscale'])

            show_legend = dplot_styles['show_legend']
            if show_legend is None:
                if (len(dplot['data']) > 1) or (np.ndim(ydata) > 1):
                    show_legend = True
            if show_legend:
                plt.legend(plot_labels, borderaxespad=0.0,
                           prop={'family': 'monospace'},
                           loc=dplot_styles['legend_loc'],
                           title=dplot_styles['legend_title'],
                           bbox_to_anchor=dplot_styles['legend_bbox'])

            if save_figures and (img_num == images_per_figure):
                if separate_figures: img_suffix = dplot_id
                else: img_suffix = str(self.fignum)

                plt.savefig(save_dir + 'image-' + img_suffix, dpi=300)

            img_num += 1

        if save_figures and (img_num < images_per_figure):
            plt.savefig(save_dir + 'image-' + str(self.fignum), dpi=300)

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

def show_results(experiment_info,
                 model=None, inv_params=None, cov=None):

    storage = DataStorage()
    data    = ResultsData()

    if model is None:
        if storage.load(experiment_info):
            data.load(storage.get('ResultsData'))
        else:
            print('      (Was computation already run?)'
                  '\nINFO: Nothing to display. Exiting.')
            exit(0)

        save_data = False
    else:
        # Disable determining of only needed calc (for inverse)
        model.measurements.run_inverse_problem_p(False)

        data.store_measurements(model.measurements)
        data.store_computation(model, model.measurements)

        if not inv_params is None: data.store_value('inv_params', inv_params)
        if not cov is None: data.store_value('cov', cov)

        from shared import get_directories
        savedir = get_directories('figs', 'mask', experiment_info)
        filename = savedir + '/' + 'results.txt'

        if not path.exists(savedir):
            makedirs(savedir)

        with open(filename, 'w') as f:
            print_status(data, filename)

        save_data = True

    dplots = DPlots(experiment_info)

    if data.store_references(dplots.get_references(), model) or save_data:
        if data.get_value('experiment_info') is None:
            data.store_value('experiment_info', experiment_info)

        storage.store('ResultsData', data.dump())
        storage.save(experiment_info)

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
