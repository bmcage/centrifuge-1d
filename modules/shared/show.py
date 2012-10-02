import matplotlib.pyplot as plt
import numpy as np
import pickle
from const import FIGS_DIR, PLOTSTYLE_ININAME
from os import makedirs, path
from shared import get_directories
from config import parse_value
try:
    import ConfigParser as configparser
except:
    import configparser


def display_table(t_measured=None, t_computed=None,
                  wl_out1_measured=None, wl_out1_computed=None,
                  gc1_measured=None, gc1_computed=None,
                  rm1_measured=None, rm1_computed=None,
                  l0_measured=None, l1_measured=None, l1_computed=None,
                  fignum=10):

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
    plt.show(block=False)

    input('Press ENTER to continue...')

def nd2strlist(nd):
    result = []
    for value in nd:
        result.append(str(value))
    return result


def has_data(x):
    if x is None:
        return False
    elif isinstance(x, np.ndarray):
        return not (x.size == 0)
    else:
        return bool(x)

dg_label_time = "Time [min]"
dg_label_length = "Sample length $L$ [cm]"
DG_AXES_LABELS = {'h': (dg_label_length, "Piezometric head $h$ [cm]"),
                  'u': (dg_label_length, "Relative saturation $u$"),
                  'MO': (dg_label_time, "Expelled water [cm]"),
                  'MI': (dg_label_time, "Inflow water [cm]"),
                  'GC': (dg_label_time, "Gravitational center [cm]"),
                  'RM': (dg_label_time, "Rotational momentum [kg.m.s$^{-1}$]"),
                  's1': (dg_label_time, "Interface s1 [cm]"),
                  's2': (dg_label_time, "Interface s2 [cm]"),
                  'WM': (dg_label_time, "Water mass [cm]"),
                  'theta': ("Water content $\\theta$", "Pressure $p$ [Pa]")}
DG_PAIRS = (('h', 'u'), ('MI', 'MO'), ('GC', 'RM'), ('s1', 's2'))

def mk_status_item(data_id, data_computed, data_measured = []):
   return {'id': data_id, 'data': (data_computed, data_measured)}

def display_status(data_plots=None):

    def compare_data(name, data_computed, data_measured, relerror, abserror):
        name_len = len(name)
        disp_all = (not data_measured is None)

        i0 = 0
        in_row = 10
        remaining = np.alen(data_computed)

        print('\n')
        while remaining > 0:
            if remaining > in_row:
                disp_items = in_row
            else:
                disp_items = remaining

            print('%s measured: ' % name,
                  disp_items * '% 10.6f' % tuple(data_measured[i0:i0+disp_items]))
            if disp_all:
                print('%s computed: ' % name,
                      disp_items * '% 10.6f' % tuple(data_computed[i0:i0+disp_items]))
                print('AbsError: ', name_len * ' ',
                      disp_items * '% 10.6f' % tuple(abs_error[i0:i0+disp_items]))
                print('Error (%):', name_len * ' ',
                      disp_items * '% 10.2f' % tuple(relerror[i0:i0+disp_items]))

            remaining = remaining - disp_items
            print(108 * '-')
            i0 = i0 + in_row

        print('LSQ error:', np.sum(np.power(data_computed - data_measured, 2)))

    if data_plots:
        for plot in data_plots:
            plot_id = plot['id']

            data_items_nr = len(plot['data'])
            if (data_items_nr == 0 ) or (not has_data(plot['data'][0])):
                continue

            if (data_items_nr == 1 ) or (not has_data(plot['data'][1])):
                value_computed = plot['data'][0]
                value_measured = rel_error = abs_error = None
            else:
                value_computed = np.asarray(plot['data'][0])
                value_measured = np.asarray(plot['data'][1])
                norm_measured = value_measured[:]
                norm_measured[norm_measured == 0.0] = 1.0e-10
                rel_error = ((value_computed - value_measured)
                             / norm_measured * 100.)

                abs_error = np.abs(value_computed - value_measured)

            compare_data(plot_id, value_computed, value_measured,
                         rel_error, abs_error)

# ResultData: hold the data of the computation
# Structure: {'lines': lines_structure, 'inv_params': inv_params, 'cov': cov}
# where:
# lines_structure: dictionary of types: {line_type: line_type_data} where:
#     line_type is (by default) one of ['computed', 'measured', 'references']
#     line_type_data is a dictionary of types: {line_id: line_data}, with
#         line_id the ID of the line
#         line_data a dict of types {data_type: (xdata, ydata)},
#              where data_type is (by default) one of
#              ['h', 'u', 'GC', 'RM', 'WM', 's1','s2'], xdata is the x-axis
#              coordinate and ydata is the y-axis coordinate (computed value)
# inv_params: dict of inverse parameters {param_name: param_value} (or None)
# cov: the covariance matrix (or None)
class ResultsData():
    def __init__(self):
        self._data = {name: None for name in ['lines', 'inv_params', 'cov']}

    def extract(self, extract_data_fn, model, referencing_parameters=[],
                 measurements=None):
        # add computed data
        (flag, value) = extract_data_fn(model)
        if not flag:
            print('Computation was not successfull. No data will be saved.')
            self._data['lines'] = None
            return

        data = {'computed': {'computed': value}}

        # add referenced data
        references = {}
        ref_num = 1
        if referencing_parameters:
            if type(referencing_parameters) == dict: # single reference
                referencing_parameters = (referencing_parameters, )

            iterable_params =  model._iterable_parameters
            for ref in referencing_parameters:
                if 'id' in ref:
                    ref_id = ref['id']
                    del ref['id']
                else:
                    ref_id = 'ref-' + str(ref_num)
                    ref_num +=1

                iters = [val for val in ref.keys() if val in iterable_params]
                if iters:
                    print('Referencing model cannot set iterable '
                          'parameters of original model:', iters)
                    exit(1)

                backup_params = model.get_parameters(ref.keys()) # backup
                model.set_parameters(ref)

                (flag, extracted_data) = extract_data_fn(model)

                model.set_parameters(backup_params) # restore

                if not flag: continue

                references[ref_id] = extracted_data

            if references:
                data['references'] = references

        # add measured data
        if measurements:
            data['measured'] = {'measured': measurements}

        self._data['lines'] = data

    def has_data(self, data_type='lines'):
        return not self._data[data_type] is None

    def get_value(self, data_name):
        if data_name in ['cov', 'inv_params', 'lines']:
            return self._data[data_name]

    def add_value(self, inv_params = None, cov = None):
        if has_data(inv_params):
            self._data['inv_params'] = inv_params
        if has_data(cov):
            self._data['cov'] = cov

    def dump(self, experiment_info):
        if not self.has_data():
            print('No data is provided. Skipping data dumping.')
            return

        savedir = get_directories('figs', 'mask', experiment_info)
        if not path.exists(savedir):
            makedirs(savedir)

        with open(savedir + 'data_results.dat', 'wb') as f:
            pickle.dump(self._data, f, pickle.HIGHEST_PROTOCOL)

    def load(self, experiment_info):
        pathdir = get_directories('figs', 'mask', experiment_info)
        filename = pathdir + 'data_results.dat'
        if not path.exists(filename):
            print('File with computation results does not exist:', filename)
            return False

        with open(filename, 'rb') as f:
            self._data = pickle.load(f)

        return True

    def get_linedatatypes(self):
         # line_type 'computed' with ID 'computed' has to be present
        return self._data['lines']['computed']['computed'].keys()
    def get_linetypes(self):
        return self._data['lines'].keys()

    def get_linedata(self, line_type, line_id):
        data = self._data['lines']
        if (line_type in data) and (line_id in data[line_type]):
            return data[line_type][line_id]
        else:
            return None

    def iterate_lines(self):
        if not self.has_data('lines'):
            yield None
        else:
            for (data_type, lines) in self._data['lines'].items():
                for (line_id, line_data) in lines.items():
                    yield (data_type, line_id, line_data)

class PlotStyles():
    def __init__(self, experiment_info):
        self._userstyles = self.load_userstyles(experiment_info)
        self._display_options = self._mk_display_options()
        self._dplotstyles = {}

    def load_userstyles(self, experiment_info):
        def get_filenames(experiment_info):
            (search_dirs, masks_dir, mask_dir) = \
              get_directories('figs', ['search', 'masks', 'mask'],
                              experiment_info)
            # include masks dir into searchable directories...
            search_dirs.append(masks_dir)
            # ...just as the single mask directory
            if experiment_info['mask']:
                search_dirs.append(mask_dir)

            filter_existing = \
              lambda fnames: list(filter(lambda fname: path.exists(fname), fnames))
            prefix_with_paths = \
              lambda fname, dirs: map(lambda cfgdir: cfgdir + fname, dirs)

            plotstyles_files = \
              filter_existing(prefix_with_paths(PLOTSTYLE_ININAME, search_dirs))

            return plotstyles_files

        def read_plotstyles_cfg(filename):
            result = {}

            # read config files
            parser   = configparser.ConfigParser()
            try:
                read_files = parser.read(filename)
            except configparser.DuplicateOptionError as E:
                print(E)
                exit(0)
            # Write data from parser to configuration
            for psection in parser.sections():
                section = psection.lower() # we ignore sections

                for option in parser.options(psection):
                    raw_value = parser.get(psection, option).strip()

                    result[option] = parse_value(raw_value)

            return result

        def deep_dictupdate(d1, d2):
            for (vkey, vvalue) in d2.items():
                if ((vkey in d1) and (type(d1[vkey]) == dict)
                    and (type(vvalue) == dict)):
                    deep_dictupdate(d1[vkey], vvalue)
                else:
                    d1[vkey] = vvalue

            return d1

        # Function body
        plotstyles_filenames = get_filenames(experiment_info)

        if not plotstyles_filenames: return {}

        plot_cfg = read_plotstyles_cfg(plotstyles_filenames[0])

        for fname in plotstyles_filenames[1:]:
            deep_dictupdate(plot_cfg, read_plotcfg(fname))

        self._userstyles = plot_cfg

        return plot_cfg

    def _get_value(self, key):
        if key in self._userstyles:
            return self._userstyles[key]
        else:
            return None

    def _mk_display_options(self):
        opts = {'separate_figures': False, 'show_figures': True,
                'save_figures': True}
        user_opts = self._get_value('options')
        if user_opts: opts.update(user_opts)

        return opts

    def get_display_options(self):
        return self._display_options

    def _mk_dplotstyles(self, dplot_id):
        dplot_styles = {tag: None for tag in ['axes_labels', 'xscale', 'yscale',
                                              'show_legend', 'legend_title',
                                              'legend_bbox', 'legend_loc',
                                              'show']}

        user_styles = self._get_value('datasets')
        if user_styles and (dplot_id in user_styles):
            dplot_styles.update(user_styles[dplot_id])

         # try to set default value
        if (dplot_styles['axes_labels'] is None) and (dplot_id in DG_AXES_LABELS):
           dplot_styles['axes_labels'] = DG_AXES_LABELS[dplot_id]

        if dplot_id in ['h', 'u']:
            if ((dplot_styles['legend_bbox'] is None)
                and (dplot_styles['legend_loc'] is None)):
                dplot_styles['legend_bbox'] = (1.01, 1.)
                dplot_styles['legend_loc'] = 2
            if dplot_styles['legend_title'] is None:
                dplot_styles['legend_title'] = 'Time [min]'

        if ((dplot_id == 'h') and (dplot_styles['show_legend'] is None)
            and (not self.get_display_options()['separate_figures'])):
            dplot_styles['show_legend'] = False

        if dplot_id == 'theta':
            if dplot_styles['yscale'] is None:
                dplot_styles['yscale'] = 'log'
            if dplot_styles['legend_loc'] is None:
                dplot_styles['legend_loc'] = 1

        if dplot_styles['legend_loc'] is None:
            dplot_styles['legend_loc'] = 4

        return dplot_styles

    def get_dplotstyles(self, dtype):
        if not dtype in self._dplotstyles:
            self._dplotstyles[dtype] = self._mk_dplotstyles(dtype)

        return self._dplotstyles[dtype]

    def _mk_linestyle(self, linetype, line_id, data_types):
        if linetype == 'measured':
            lineopt = 'x'
        else:
            lineopt = '.'

        line_styles = {dtype: lineopt for dtype in data_types}

        if 'h' in data_types: line_styles['h'] = '-'
        if 'u' in data_types: line_styles['u'] = '-'
        if ('theta' in data_types) and (not linetype == 'measured'):
             line_styles['theta'] = '-'

        user_styles = self._get_value('lines')
        if user_styles and (line_id in user_styles):
            line_styles.update(user_styles[line_id])

        return line_styles

    def get_linestyle(self, line_type, line_id, data_types):
        return self._mk_linestyle(line_type, line_id, data_types)

class DPlots():
    def __init__(self, data, experiment_info):
        self._data = data
        self._experiment_info = experiment_info
        self._dplots = None
        self._plotstyles = None

    def display(self, fignum = 1):

        def _mk_dplots_bucket(data_types, plot_styles):
            dplots_bucket = \
              {dtype: {'id': dtype, 'data': [],
                       'styles': plot_styles.get_dplotstyles(dtype)}
                for dtype in data_types}

            return dplots_bucket

        def _add_plotline(line_type, line_id, line_data, data_types,
                          plot_styles, dplots_bucket):

            line_styles = plot_styles.get_linestyle(line_type, line_id,
                                                        data_types)
            if 'label' in line_styles:
                label = line_styles['label']
            else:
                label = line_id

            for (data_type, data_value) in line_data.items():
                if not data_value: continue

                # we skip other 'h' and 'u' data, as it would be mess
                if (data_type in ['h', 'u']) and (not label == 'computed'):
                    continue

                (xdata, ydata) = (data_value[0], data_value[1])

                if not (has_data(xdata) and has_data(ydata)): continue

                if ((data_type in ['h', 'u']) and (len(data_value) > 2)
                    and has_data(data_value[2]) and (not 'label' in line_styles)):
                    ilabel = ['% 6d' % (ti/60.) for ti in data_value[2]]
                else:
                    ilabel = label

                item = (xdata, ydata, ilabel, line_styles[data_type])

                dplots_bucket[data_type]['data'].append(item)

        def _filter_dplots(dplots_bucket):
            for name in list(dplots_bucket.keys()):
                if dplots_bucket[name]['styles']['show'] == False:
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

        def _mk_dplots(data, experiment_info):
            plot_styles   = PlotStyles(experiment_info)
            data_types    = data.get_linedatatypes()
            dplots_bucket = _mk_dplots_bucket(data_types, plot_styles)

            for (line_type, line_id, line_data) in data.iterate_lines():
                _add_plotline(line_type, line_id, line_data, data_types,
                              plot_styles, dplots_bucket)

            ordered_dplots = _order_dplots(_filter_dplots(dplots_bucket))

            return (ordered_dplots, plot_styles)

        def _show_dplots(ordered_dplots, display_options, experiment_info):
            nonlocal fignum

            separate_figures = display_options['separate_figures']
            save_figures     = display_options['save_figures']

            if save_figures:
                experiment_info = experiment_info
                if experiment_info['mask']:
                    figs_dir_type = 'mask'
                else:
                    figs_dir_type = 'data'

                save_dir = get_directories('figs', figs_dir_type,
                                           experiment_info)

                if not path.exists(save_dir):
                    makedirs(save_dir)


            if separate_figures:
                images_per_figure = 1
            else:
                images_per_figure = 6

            fignum -= 1
            img_num = 2^20 # high initialization, so that first fig is created

            for dplot in ordered_dplots:
                if not dplot['data']: continue
                dplot_id = dplot['id']

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
                    plt.subplot(3,2,img_num)

                # plot the supplied data
                plot_labels = []
                for dplot_data in dplot['data']:
                    (xdata, ydata, line_label, plot_style) = dplot_data
                    plt.plot(xdata, ydata, plot_style)
                    if type(line_label) == str:
                        plot_labels.append(line_label)
                    else:
                        plot_labels.extend(line_label)

                dplot_styles = dplot['styles']
                (xlabel, ylabel) = dplot_styles['axes_labels']
                plt.xlabel(xlabel)
                plt.ylabel(ylabel)
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
                    else: img_suffix = str(fignum)

                    plt.savefig(save_dir + 'image-' + img_suffix, dpi=300)

                img_num += 1

            if save_figures and (img_num < images_per_figure):
                plt.savefig(save_dir + 'image-' + str(fignum), dpi=300)

        def _show_status(data):
            status_items = []

            measurements = data.get_linedata('measured', 'measured')
            computed     = data.get_linedata('computed', 'computed')

            for (key, m_data) in measurements.items():
                if m_data[1] is None: continue

                if key in computed:
                    if key == 'theta':
                        c_value = computed[key][2]
                        m_value = m_data[0]
                    else:
                        c_value = computed[key][1][1:]
                        m_value = m_data[1]

                if c_value is not None:
                    status_items.append(mk_status_item(key, c_value, m_value))

            if status_items:
                display_status(status_items)

        # function body
        data   = self._data

        if self._dplots is None: # if data is provided, generate dplots
            if not data.has_data('lines'):
                print('No data is provided. Nothing to display.')
            else:
                (self._dplots, self._plotstyles) = \
                  _mk_dplots(data, self._experiment_info)

        if not self._dplots is None: # is True if data was provided
            _show_status(data)

            if data.has_data('cov'):
                print('\nCov:\n', data.get_value('cov'))
            if  data.has_data('inv_params'):
                params = data.get_value('inv_params')
                print('\nOptimal parameters found:')
                for (name, value) in params.items():
                    if name == 'ks':
                        print('  Ks [cm/s]: {: .8g}'.format(value))
                    else:
                        print('  {:9}: {: .8g}'.format(name, value))

            _show_dplots(self._dplots, self._plotstyles.get_display_options(),
                         self._experiment_info)
            plt.show(block=False)

        input('Press ENTER to continue...')
