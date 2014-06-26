from __future__ import print_function, division

import numpy as np, matplotlib, pickle
import sys
from ...const import PLOTSTYLE_ININAME, DUMP_DATA_VERSION, DUMP_DATA_FILENAME
from os import makedirs, path
from ...shared import get_directories, parse_value, yn_prompt, filter_indices
from ...config import ModulesManager, load_model
from collections import OrderedDict, defaultdict
from .functions import has_data, compare_data
from .saturation_curve import create_SC

if sys.version_info[0] == 2:
    #python 2.7
    import ConfigParser as configparser
    input = raw_input
else:
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

# Linestyle is searched in order from most specific to least specific and return
# the first found value of specified option. The order is:
#   1. 'line_id'   -> 'fig_id'
#   2. 'line_id'   -> '_base_'
#   3. '_default_' -> 'fig_id'
#   4. '_default_' -> '_base_'
# Linestyle options:
#    Assigned only to '_base_' of line_id:
#        'xdata', 'ydata', 'label', 'legend_data', 'order',
#    Assigned anywhere:
#        'width', 'symbolsize', 'lineopt'
LINESTYLES_DEFAULT = {'_default_': {'_base_': {'width': 1, 'symbolsize': 6,
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
#    'show', 'show_legend', 'legend_title', 'legend_bbox', 'legend_loc',
#    'overlay_alpha', 'overlay_color',
#    'overlay_x' = None or (x0, x1), overlay_y = None or (y0, y1))
#        - if both are None, no overlay is displayed
#        - x0, x1, y0, y1 are float or None (= use default value)

FIG_OPTIONS_DEFAULTS = {'show': True, 'legend_loc': 4, 'order': 999,
                        'show_legend': True, 'overlay_color': 'blue',
                        'overlay_x': None, 'overlay_y': None,
                        'overlay_alpha': 0.2}

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
                      'xlabel': dg_label_time,
                      'ylabel': "Inflow water [{}]",
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
                      'xlabel': dg_label_time,
                      'ylabel': "Water mass balance [{}]",
                      'xtype':  'time',
                      'ytype':  'length'},
               'WM_in_tube':  {'title':  'Water mass in tube',
                               'xlabel': dg_label_time,
                               'ylabel': "Water mass in tube [{}]",
                               'xtype':  'time',
                               'ytype':  'length'},
               'theta':   {'title': 'Pressure $p$ vs $\\theta$',
                           'xlabel': "Water content $\\theta${}",
                           'ylabel': "Pressure $p$ [{}]",
                           'xtype':  'none',
                           'ytype':  'pressure',
                           'yscale': 'log',
                           'legend_loc': 1,
                           'overlay_x': (None, None)},
               'relsat':  {'title':  'Negative hydraulic head $h$',
                           'xlabel': "Effective saturation $S_e$ [{}]",
                           'ylabel': "Negative hydraulic head $h$ [{}]",
                           'xtype':  'none',
                           'ytype':  'length',
                           'yscale': 'log',
                           'legend_loc': 1,
                           'overlay_x': (None, None)},
               'K':       {'title':  'Hydraulic conductivity $K(\\theta)$',
                           'xlabel': "Water content $\\theta${}",
                           'ylabel': "Hydraulic conductivity $K(\\theta)$ [{}]",
                           'xtype':  'none',
                           'ytype':  'velocity',
                           'yscale': 'log',
                           'overlay_x': (None, None)},
               'K_u':     {'title':  'Hydraulic conductivity $K(S_e)$',
                           'xlabel': "Effective saturation $S_e${}",
                           'ylabel': "Hydraulic conductivity $K(S_e)$ [{}]",
                           'xtype':  'none',
                           'ytype':  'velocity',
                           'yscale': 'log',
                           'overlay_x': (None, None)},
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

FIGURES_PAIRS = (('h', 'u'), ('MI', 'MO'), ('GC', 'RM'), ('K', 'K_u'),
                 ('gF_MT', 'gF_MO'), ('dgF_MT', 'dgF_MO'), ('s1', 's2'),
                 ('theta', 'relsat'))

FIGURES_IDS = list(FIGURES_DEFAULTS.keys())

DISPLAY_OPTIONS = {'separate_figures': False, 'show_figures': True,
                   'save_figures': True, 'matplotlib_backend': None,
                   'show_figures_titles': None, 'comparison_table': False,
                   'save_formats': ['png'], 'figures_dpi': 92}

def set_default_units(figures_styles):
    for fig_style in figures_styles.values():
        fig_style['xunit'] = DEFAULT_UNITS[fig_style['xtype']]
        fig_style['yunit'] = DEFAULT_UNITS[fig_style['ytype']]

    return figures_styles

def get_unit_coef(unit_base):
    """
      Return the coeficient for the new unit type to be used, so that
      the internal unit is converted to the new unit.
    """
    default_units =  ['cm', 's', 'pa', 'g', 'gf', 'cm/s', 'rpm', '']

    unit = unit_base.lower()
    # units used for computation are: cm, s, pa, gf and "no units"
    if unit in default_units: coef = 1.0
    elif unit == 'mm': coef = 10.
    elif unit in 'min': coef = 1./60.
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
    line_ids = (line_id, '_default_')
    fig_ids = (fig_id, '_base_') if fig_id is not None else ('_base_', )

    for line_name in line_ids:
        for fig_name in fig_ids:
            if ((line_name in linestyles)
                 and (fig_name  in linestyles[line_name])
                 and (name in linestyles[line_name][fig_name])):
                return linestyles[line_name][fig_name][name]

    return not_found

def get_figure_option(figstyles, fig_id, name, not_found=None):
    fig_style = figstyles[fig_id]

    if name in fig_style:
        value = fig_style[name]
    elif name in FIG_OPTIONS_DEFAULTS:
        value = FIG_OPTIONS_DEFAULTS[name]
    else:
        value = not_found

    return value

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

    if (not computed) or (not measurements):
        return

    stream = None if filename is None else open(filename, 'w')

    try:
        total_LSQ_error = 0.0
        # compare measured vs. computed data
        for (key, m_data) in measurements.items():
            if (m_data[1] is None) or (not key in computed):
                continue

            if key == 'theta':
                if  len(computed[key]) == 2:
                    # OEPS TODO, fit_retention_curve, this is working
                    #   when writing to file !!! WHY ??
                    continue
                else:
                    c_value = computed[key][3]
            else:
                c_value = computed[key][1]
            m_value = m_data[1]

            if c_value is not None:
                LSQ_error, RMS_error = compare_data(key, c_value, m_value, stream)
                total_LSQ_error += LSQ_error

        print('\nTotal LSQ error: ', total_LSQ_error, file=stream)

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
                    if np.iterable(value):
                        print('  Ks [cm/s]: {0!s}'.format(value), file=stream)
                    else:
                        print('  Ks [cm/s]: {: .8g}'.format(value), file=stream)
                elif np.iterable(value):
                    print('  {0:9}: {1!s}'.format(name, value), file=stream)
                else:
                    print('  {:9}: {: .8g}'.format(name, value), file=stream)

        params = data.get('print_params')
        if not params is None:
            print('\nUsed parameters:', file=stream)
            for name in params:
                value = data.get(name)
                print (name, value)
                if value is None:
                    continue
                if name == 'ks':
                    if np.iterable(value):
                        print('  Ks [cm/s]: {0!s}'.format(value), file=stream)
                    else:
                        print('  Ks [cm/s]: {: .8g}'.format(value), file=stream)
                elif np.iterable(value):
                    print('  {0:9}: {1!s}'.format(name, value), file=stream)
                else:
                    print('  {:9}: {: .8g}'.format(name, value), file=stream)

    finally:
        if not filename is None: stream.close()


    #output synthetic data for direct runs
    streamsynthetic = None if filename is None else \
                            open(filename[:-5]+"_synth.ini", 'w')
    if not (streamsynthetic is None):
        #section heading
        print('[experiment]\n', file=streamsynthetic)

    try:
        for (key, c_data) in computed.items():
            if key == 'theta':
                if  len(computed[key]) == 2:
                    # OEPS TODO, fit_retention_curve, this is working
                    #   when writing to file !!! WHY ??
                    continue
                else:
                    c_value = computed[key][3]
            else:
                c_value = computed[key][1]
            c_xvalue = c_data[0]

            if c_value is not None and key in data._styles['output_computed_as_meas']:
                output_synthetic(key, c_xvalue, c_value, streamsynthetic)
    finally:
        if not filename is None: streamsynthetic.close()

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

def is_dict_p(value):
    return type(value) in (dict, OrderedDict, defaultdict)

def deep_dictupdate(d1, d2):
    """ Recursively update dictionary d1 with values of dictionary d2. """

    for (vkey, vvalue) in d2.items():
        if (vkey in d1) and is_dict_p(d1[vkey]) and is_dict_p(vvalue):
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
       any data. Transform 'overlay_x' and 'overlay_y' if needed.
     """
    # Resolve 'overlay_x' and 'overlay_y'
    figures_styles = styles['figures']
    line_data = data.get_linedata('computed')
    if 'u' in line_data:
        udata = line_data['u'][1][:,:].flatten()
        udata = udata[udata != 0.]
        u_min = np.min(udata)
        u_max = np.max(udata)
    else:
        #output has no u, eg saturated data!
        u_min = 0
        u_max = 1
    for fig_id in displayed_figs:
        overlay_x = get_figure_option(figures_styles, fig_id, 'overlay_x')
        overlay_y = get_figure_option(figures_styles, fig_id, 'overlay_y')

        if (overlay_x is None) and (overlay_y is None):
            continue

        (ox0, ox1) = (None, None) if overlay_x is None else overlay_x
        (oy0, oy1) = (None, None) if overlay_y is None else overlay_y

        if fig_id in ('K_u', 'relsat'):
            ox0 = u_min if ox0 is None else ox0
            ox1 = u_max if ox1 is None else ox1
        elif fig_id in ('K', 'theta'):
            if 'theta' in line_data:
                theta_s = line_data['theta'][1][0] # = theta at full saturation
                # normally theta_r << 1 hence it won't be visible in the figure
                # if we consider theta_r == 0
                ox0 = theta_s*u_min if ox0 is None else ox0
                ox1 = theta_s*u_max if ox1 is None else ox1

        figures_styles[fig_id]['overlay_x'] = (ox0, ox1)
        figures_styles[fig_id]['overlay_y'] = (oy0, oy1)

    # build a list of figs that actually contain some data
    nonempty_figs = []

    if not displayed_figs:
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
            if ((not fig_id in line_data)
                 or (fig_id in ['h', 'u']) and (not line_id == 'computed')):
                continue

            line_value = line_data[fig_id]

            # Process also special 'legend_data' (e.g. for 'h' and 'u' it's times)
            if len(line_value) > 2:
                (xdata, ydata, legend_data) = line_value
                if type(legend_data) is np.ndarray:
                    legend_data = ['% 6.1f' % (ti/60.) for ti in legend_data]
            else:
                (xdata, ydata) = line_value
                legend_data = None

            if (not has_data(xdata)) or (not has_data(ydata)): continue

            # This fig surely will have some data
            if not fig_id in nonempty_figs:
                nonempty_figs.append(fig_id)

            if not fig_id in line_style:
                line_style[fig_id] = {}

            line_fig_style = line_style[fig_id]

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
                    filter_idxs = np.zeros((filter_size, ), dtype=bool)

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
                # switch xdata <-> ydata
                (xdata, ydata) = (ydata, xdata)

            line_fig_style['xdata'] = xdata
            line_fig_style['ydata'] = ydata

    return nonempty_figs

def get_shown_figs_ids(figures_styles):
    # Keep only figures with 'show' == True
    figs_ids = [fig_id for fig_id in figures_styles.keys()
                  if get_figure_option(figures_styles, fig_id, 'show')]

    return figs_ids

def draw_figure(fig, fig_id, figs_styles, lines_ids, lines_styles,
                show_titles):

    import matplotlib.pyplot as plt
     # plot the lines
    legend_labels = []

    xunit = get_figure_option(figs_styles, fig_id, 'xunit')
    yunit = get_figure_option(figs_styles, fig_id, 'yunit')

    for line_id in lines_ids:
        xdata = get_line_option(lines_styles, line_id, 'xdata', fig_id)
        ydata = get_line_option(lines_styles, line_id, 'ydata', fig_id)

        if (xdata is None) or (ydata is None):
            continue

        width = get_line_option(lines_styles, line_id, 'width', fig_id)
        symbolsize = get_line_option(lines_styles, line_id,
                                             'symbolsize', fig_id)
        xcoef = get_unit_coef(xunit)
        ycoef = get_unit_coef(yunit)
        if not xcoef == 1.0:
            xdata = xcoef * np.asarray(xdata, dtype=float)
        if not ycoef == 1.0:
            ydata = ycoef * np.asarray(ydata, dtype=float)

        plot_style = get_line_option(lines_styles, line_id, 'lineopt',
                                             fig_id)
        if (len(xdata.shape) > 1) and np.iterable(plot_style):
            max_styles = len(plot_style)
            for ind in range(len(xdata[0])):
                entryx = xdata[:, ind]
                entryy = ydata[:, ind]

                plt.plot(entryx, entryy, plot_style[ind%max_styles],
                         linewidth=width, markersize=symbolsize)
        else:
            plt.plot(xdata, ydata, plot_style, linewidth=width,
                     markersize=symbolsize)

        # Extend the legend labels
        legend_label = get_line_option(lines_styles, line_id,
                                       'legend_data', fig_id)
        if legend_label is None:
            legend_label = get_line_option(lines_styles, line_id,
                                           'label', fig_id, line_id)
        if np.isscalar(legend_label):
            legend_labels.append(legend_label)
        elif legend_label is not None:
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

    #handle overlays
    overlay_x = get_figure_option(figs_styles, fig_id, 'overlay_x')
    overlay_y = get_figure_option(figs_styles, fig_id, 'overlay_y')
    if (not overlay_x is None): # either both are None or none
        (ox0, ox1) = (overlay_x[0] or xmin or 0.0, overlay_x[1] or xmax or 1)
        (oy0, oy1) = (overlay_y[0] or ymin or 0.0, overlay_y[1] or ymax or 1e10)

        ax = plt.gca()
        ax.fill_between(np.array([ox0, ox1],float), oy0, oy1,
                        facecolor=get_figure_option(figs_styles, fig_id,
                                                    'overlay_color'),
                        alpha=get_figure_option(figs_styles, fig_id,
                                                 'overlay_alpha'))

    if show_titles:
        plt.suptitle(get_figure_option(figs_styles, fig_id, 'title'))

    show_legend = get_figure_option(figs_styles, fig_id, 'show_legend')
    if show_legend:
        legend_loc   = get_figure_option(figs_styles, fig_id, 'legend_loc')
        legend_title = get_figure_option(figs_styles, fig_id, 'legend_title')
        legend_bbox  = get_figure_option(figs_styles, fig_id, 'legend_bbox')
        plt.legend(legend_labels, borderaxespad=0.0, title=legend_title,
                   prop={'family': 'monospace'}, loc=legend_loc,
                   bbox_to_anchor=legend_bbox)


################################################################################
#                                Data storage                                  #
################################################################################

class DataStorage:
    """
      Object for holding computed, measured data and additional data
      and displays them applying user styles.
    """

    def __init__(self, experiment_info):
        self._experiment_info = experiment_info
        self._modman = None
        self._data = {'lines': {}}                   # stored data

        # Default values (merged with plotstyles inifiles)
        styles = {'options': DISPLAY_OPTIONS, 'lines': LINESTYLES_DEFAULT,
                  'figures': set_default_units(FIGURES_DEFAULTS),
                  'plots_keep': {}, 'plots_remove': {}, 'params_ref': {},
                  'lines_order': (),
                  'output_computed_as_meas': {},
                 }

        # Read user plotysles inifiles
        plotstyles_filenames = get_filenames(experiment_info)

        user_styles = {}
        if plotstyles_filenames:
            for fname in plotstyles_filenames:
                deep_dictupdate(user_styles, read_plotstyles_file(fname))

        update_styles(styles, user_styles)

        self._styles = styles

        display_options = styles['options']

        if 'matplotlib_backend' in display_options:
            matplotlib_backend = display_options['matplotlib_backend']
        if matplotlib_backend: # chosen other backend than the default one
            matplotlib.use(matplotlib_backend)

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

#        print("Storing Computation %s with main parameters:" % ID,
#                  model.get_parameters(['ks', 'hi', 'ki', 'ui', 'n', 'gamma']),
#                    ', h_init=', getattr(model, 'h_init', None))
        if self._modman is None:
            self._modman = ModulesManager()

        solver_module = self._modman.find_module(model.exp_type,
                                                 submodule='run')
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
            if ID == 'computed':
                #no sense continuing if computed does not work
                raise Exception('Could not obtain computed data')
        else:
            line_style = self._styles['lines'][ID]
            data = {}
            if self.get('print_params') is None:
                #store from model extra param to print from computation
                self.store('print_params', model.get_value('print_params'))
            if self.get('print_params')  and ID == 'computed':
                for name in self.get('print_params'):
                    val = model.get_parameters([name])
                    if name in val:
                        val = val[name]
                    else:
                        val = None
                    if val in [None, '']:
                        val = getattr(model, name, None)
                    self._data[name] = val

            # Store all computed data
            for (name, xvalue, yvalue) in measurements.iterate_calc_measurements():
                if xvalue is None or yvalue is None:
                    continue
                # make a local copy as array may be overwritten
                xvalue = xvalue.copy()
                yvalue = yvalue.copy()

                if name in ('h', 'u'):
                    t = measurements.get_times()
                    xvalue = xvalue.transpose()
                    yvalue = yvalue.transpose()
                    data[name] = (xvalue, yvalue, t)
                else:
                    data[name] = (xvalue, yvalue)

                line_style['xdata'] = xvalue
                line_style['ydata'] = yvalue

            # Store extra data
            # a) Retention curve based on theta
            if hasattr(model, 'SC'):

                SC = model.SC
                if SC.canrefine_h() and ID == 'computed':
                    self._data['hi'] = SC.refinable_h()

                theta_s = (model.theta_s if hasattr(model, 'theta_s')
                            else  model.porosity)
                theta_r = model.theta_r if hasattr(model, 'theta_r') else 0.0

                (p, theta) = SC.retention_curve(theta_s, model.density, model.g,
                                                theta_r=theta_r)
                (rc_h, rc_u) = SC.retention_curve_u(g=model.g,
                                                    rho=model.density)

                if hasattr(model, 'p') or hasattr(model, 'pc'):
                    # TODO: can we not remove this?? Is stored in measured NOW!
                    # returning more than lenth two in data no longer supported !!
                    p_meas = np.asarray(model.p if hasattr(model, 'p')
                                        else model.pc, dtype=float)

                    (p_user, theta_user) = \
                      SC.retention_curve(theta_s, model.density, model.g,
                                         theta_r=theta_r, p=p_meas)
                    (rc_h_user, rc_u_user) = \
                      SC.retention_curve_u(g=model.g, rho=model.density,
                                           p=p_meas)

                    data['theta']  = (p, theta)
                    #, p_user, theta_user)   # no longer supported
                    data['relsat'] = (-rc_h, rc_u)
                    #, -rc_h_user, rc_u_user)   # no longer supported
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
                        data['K'] = thK
                        data['K_u']  = uK
                    except:
                        import traceback
                        print(traceback.format_exc())

            self._data['lines'][ID] = data
            self.store('experiment_info', model.experiment_info)

        return flag

    def store_references(self, model=None):
        """ Store computations corresponding to referenced parameters. """

        stored_references = self.get('references', not_found={})
        user_references = self._styles['params_ref']

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

        iterable_params = model._iterable_parameters
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
                # TODO: load plotstyles.ini with cfg loader to have
                #       all parameters names in lowercase
                backup_SC = model.SC
                backup_typeSC = model.SC.typeSC()
                if 'SC_type' in ref_params and not ('sc_type' in ref_params):
                    ref_params['sc_type'] = ref_params['SC_type']
                if (not ('sc_type' in ref_params)):
                    print('Referencing model does not contain "SC_type", '
                          'cannot set parameters of the saturation curve.'
                          '\nSkipping...')
                try:
                    model.SC = create_SC(ref_params)
                except Exception as E:
                    print(E)
                    print('Could not process ref_params. Continuing...')

            model.set_parameters(ref_params)

            flag = self.store_computation(model, model.measurements, ID=ref_id)

            if not flag:
                print('Reference parameters: ', ref_params)

             # restore original SC
            if hasattr(model, 'SC'):
                if backup_typeSC != model.SC.typeSC():
                    model.SC = backup_SC
            model.set_parameters(backup_params)

        self.store('references', stored_references)

        return True

    def store_measurements(self, measurements, model=None):
        """ Store measured (supplied) data. """

        for mtype in ('original', 'measured'):
            m = {name: (xvalue, yvalue) for (name, xvalue, yvalue)
                 in measurements.iterate_meas_measurements(mtype == 'original',
                                                           model)}

            self._data['lines'][mtype] = m

    def get_linedata(self, line_id, not_found=None):
        """ Get data for specified line. """
        return self._data['lines'].get(line_id, not_found)

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

        pathdir = get_directories('figs', 'mask', self._experiment_info)
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

    def display(self, fignum=1, ext='', noblock=False):
        """ Display the figures and/or write them to files. """

        import matplotlib.pyplot as plt
        styles = self._styles

        display_options = styles['options']

        show_figures     = display_options['show_figures']
        save_figures     = display_options['save_figures']
        separate_figures = display_options['separate_figures']
        show_titles      = display_options['show_figures_titles']

        if not (save_figures or show_figures):
            return

        # make sure user did not pass True
        show_titles = separate_figures if (show_titles is None) else False

        figs_styles = styles['figures']
        figs_ids = assign_data(styles, get_shown_figs_ids(figs_styles), self)

        if not figs_ids:
            return

        figs_ids = order_figures(figs_styles, figs_ids)

        lines_styles = styles['lines']
        lines_ids = styles['lines_order']

        print_status(self)

        print('\nGenerating figures... ', end='')

        if save_figures:
            experiment_info = self._experiment_info
            figs_dir_type = 'mask' if experiment_info['mask'] else 'data'

            save_dir = get_directories('figs', figs_dir_type, experiment_info)

            if not path.exists(save_dir):
                makedirs(save_dir)

            save_figs_list = []

        images_per_figure = 1 if separate_figures else 6

        fignum -= 1
        img_num = 2^20 # high initialization, so that first fig is created
        for fig_id in figs_ids:
            # resolve figure and subplot
            if img_num > images_per_figure:
                img_num = 1
                fignum += 1

                if separate_figures:
                    fig = plt.figure(fignum)
                    img_suffix = fig_id
                else:
                    fig = plt.figure(fignum, figsize=(16, 8.5))
                    plt.subplots_adjust(wspace=0.15, left=0.06, right=0.85)
                    img_suffix = str(fignum)

                if save_figures:
                    save_figs_list.append((fig, img_suffix))

            if not separate_figures:
                plt.subplot(3, 2, img_num)

            draw_figure(fig, fig_id, figs_styles, lines_ids, lines_styles,
                        show_titles)

            img_num += 1

        if display_options['comparison_table']:
            display_table(self, fignum+1)

        if show_figures:
            print('Displaying figures... ', end='')
            try:
                plt.show(block=False)
            except: # Older matplotlib compatibility
                plt.ion()
                plt.show()

        if save_figures:
            print('Saving figures... ', end='')

            save_formats = display_options['save_formats']
            if np.isscalar(save_formats):
                save_formats = (save_formats, )
            figures_dpi = display_options['figures_dpi']

            for (fig, img_suffix) in save_figs_list:
                for save_format in save_formats:
                    try:
                        fig.savefig(save_dir + 'image-' + img_suffix + ext
                                     + '.' + save_format,
                                    format=save_format, dpi=figures_dpi)
                    except:
                        print('Saving figures to format ' + save_format
                               + ' was not successful. Skipping...')

        print('Done.')              # Generating/displaying/saving figures

        if noblock:
            if show_figures:
                plt.close('all')
        else:
            if show_figures:
                input('\nPress ENTER to continue...')
                plt.close('all')

################################################################################
#                         Module's provided functions                          #
################################################################################
def show_inbetween_result(experiment_info, model, inv_params, cov, ext):
    """
    store inbetween resuls, as files with extra ext ending.
    """
    data = DataStorage(experiment_info)
    model.measurements.run_inverse_problem_p(False)

    data.store_measurements(model.measurements, model)
    data.store_computation(model, model.measurements)
    data.store('inv_params', inv_params)
    data.store('cov', cov)
    savedir = get_directories('figs', 'mask', experiment_info)
    filename = savedir + '/' + 'results' + ext + '.txt'

    if not path.exists(savedir):
        makedirs(savedir)


    print_status(data, filename)

    data.display(ext=ext, noblock=True)

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

    if data.store_references(model) or save_data:
        data.save()

    data.display()

def display_table(data, fignum=None):
    """ Display comparison table of measured vs. computed results. """

    from matplotlib.ticker import NullLocator
    import matplotlib.pyplot as plt

    majorLocator = NullLocator()

    AXESPAD = 0.02  # default (unexported) value in matplotlib.table
    NAMECELL_WIDTH = 0.1 # width of first col (containing name)

    fig = plt.figure(fignum)
    ax = fig.add_axes((0,0,1,1))
    ax.xaxis.set_major_locator(majorLocator)
    ax.yaxis.set_major_locator(majorLocator)

    colsnum = 10

    cols_widths = ([NAMECELL_WIDTH]
                   + [(1-NAMECELL_WIDTH-2*AXESPAD)/colsnum] *colsnum)

    computed     = data.get_linedata('computed')
    measurements = data.get_linedata('measured')

    rows = []
    for (key, m_data) in measurements.items():
        if (m_data[1] is None) or (not key in computed):
            continue


        c_value = computed[key][3] if key == 'theta' else computed[key][1]
        m_value = m_data[1]

        if c_value is None:
            continue

        compare_results(key, c_value, m_value, colsnum, rows.append)

    plt.table(cellText=rows, colLabels=None, rowLabels=None, lod=False,
              colWidths=cols_widths, loc=2);

def compare_results(name, c_value, m_value, colsnum, apply_row_fn,
                    padding=True):
    """
        Perform comparison of two arrays - 'c_value' and 'm_value'. Apply the
        'apply_row_fn' on one row of the comparison. Number of rows is
        determined by the length of the input arrays and number of columns
        'colsnum'. If 'padding' is True empty strings are appended to fill data
        to 'colsnum'.
        Rows 1-4 are repeated untill the end of c_value and m_value is reached:
        1.row: [name, colsnum * measured_data]
        2.row: [name, colsnum * computed_data]
        3.row: [AbsError, colsnum * abserror_data]
        4.row: [RelError, colsnum * relerror_data]

        Lastly, row with LSQ error is added (and empty row afterwards)
    """

    error = c_value - m_value
    abs_error = np.abs(error)
    m_norm = np.max(np.abs(m_value), 1e-60)
    rel_error = abs_error / m_norm * 100
    rel_error[rel_error > 1e50] = np.inf
    rel_error = np.sign(error) * rel_error

    float_display_length = 12
    fstr = '{: ' + str(float_display_length) + '.6f}' # '{ 12.6f}}'

    rows_names = ('{} measured'.format(name), '{} computed'.format(name),
                  'AbsError: ', 'Error (%):')
    rows_values = (m_value, c_value, abs_error, rel_error)

    items_count = np.alen(c_value)
    i0 = 0
    while i0 < items_count:
        i1 = min(i0 + colsnum, items_count)
        for (row_name, row_value) in zip(rows_names, rows_values):
            apply_row_fn([row_name] +
                         [fstr.format(val) for val in row_value[i0:i1]]
                         + ['']*(i0+colsnum-i1)*bool(padding))

        apply_row_fn(['LSQ error ', "'" + name + "': ",
                      fstr.format(np.sum(np.power(error, 2)))]
                       + ['']*(colsnum - 2)*bool(padding))
        apply_row_fn(['']*(colsnum + 1)*bool(padding))
        i0 = i1

def output_synthetic(name, xvalue_computed, value_computed = None,
                 stream=None):
    if stream is None: return

    print('%s =  [' % (name+'_xvalues'), file=stream, end="")
    for val in xvalue_computed[:-1]:
        print(val, file=stream, end="")
        print(', ', file=stream, end="")
    print(xvalue_computed[-1], "]", file=stream)

    print('%s =  [' % name, file=stream, end="")
    for val in value_computed[:-1]:
        print(val, file=stream, end="")
        print(', ', file=stream, end="")
    print(value_computed[-1], "]", file=stream)
