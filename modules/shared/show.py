import matplotlib.pyplot as plt
import numpy as np
from const import FIGS_DIR
from os import makedirs, path

def y2x(y, s1, s2):
    s1_len = np.alen(s1)
    if s1_len != np.alen(s2):
        print('Interfaces array ''s1'' and ''s2'' have to be of the same'
              'lenght. Cannot proceed.')
        exit(1)
    x = np.empty([s1_len, len(y)], float)

    ds = s2 - s1

    for i in range(s1_len):
        x[i, :] = s1[i] + y * ds[i]

    return x

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
                  'RC': ("Water content $\\theta$", "Pressure $p$ [Pa]")}
DG_PAIRS = (('h', 'u'), ('MI', 'MO'), ('GC', 'RM'), ('s1', 's2'))

def add_plotline(dplot, xdata, ydata, label=None, line_opts=None):
    if line_opts is None:
        item = (xdata, ydata, label)
    else:
        item = (xdata, ydata, label, line_opts)

    plot['data'].append(item)

    return item

def make_dplot(plot_id, plot_items, legend_loc=None, show_legend=None,
            legend_title=None, legend_bbox=None, xscale=None, yscale=None,
            axes_labels=None):

    plot = {'id': plot_id, 'data': plot_items}
    if not legend_loc is None:
        plot['legend_loc'] = legend_loc
    if not axes_labels is None:
        plot['axes_labels'] = axes_labels
    if not legend_title is None:
        plot['legend_title'] = legend_title
    if not legend_bbox is None:
        plot['legend_bbox'] = legend_bbox
    if not show_legend is None:
        plot['show_legend'] = show_legend
    plot['xscale'] = xscale
    plot['yscale'] = yscale

    return plot

def draw_graphs(times=None, t_ref = None, y = None, h = None, u = None,
                s1 = None, s1_ref=None, s2 = None, s2_ref=None,
                mass_out = None, mass_out_ref = None,
                mass_in = None, mass_in_ref = None,
                GC = None, GC_ref = None,
                RM = None,  RM_ref = None, WM = None, WM_ref = None,
                fignum = 1, save_figures=False, separate_figures=False,
                save_as_text=False, draw_equilibrium=False,
                show_figures=False, experiment_info=None,
                plots = None,
                model=None):

    def add_legend(lines, legend_data=None, legend_title=None, legend_loc=1,
                   legend_type='figlegend'):
        if legend_type == 'figlegend':
            legendfn = plt.figlegend
        elif legend_type == 'legend':
            legendfn = plt.legend
        else:
            raise ValueError('Unknown legend type: ', legend_type)

        legendfn(lines, legend_data, legend_loc, borderaxespad=0.0,
                 title=legend_title, prop={'family': 'monospace'})

    print('\n', 30*'-', '\n  Displaying results...\n', 30*'-')

    if not plots:
        print('No data suplied. Nothing to display.')
        exit(0)

    def narrow_plots(plots):
        conformed_plots = {key: [] for key in DG_AXES_LABELS.keys()}

        # item: {'id':string, 'data':item_data_list
        #        [, 'labels':(xlabel, ylabel)][, 'legend_loc':loc]
        #        [, 'legend_bbox': bbox_desc][, 'legend_title': string]
        #        [, 'show_legend': bool]}
        # item_data: (x, y[, label[, draw_opts]])
        for item in plots:
            item_id   = item['id']
            item_data = item['data']

            if not (item_data and has_data(item_data[0][1])): continue

            # create a duplicate of existing values...
            new_item = {}
            for (key, value) in item.items():
                if key in ['data']: continue
                new_item[key] = value

            #  ..and fill optional values:

            if not 'label' in item:            # fill axes labels
                new_item['axes_labels'] = DG_AXES_LABELS[item_id]


            if not 'legend_loc' in item:       # fill legend location
                if item_id in ['h', 'u']:
                    new_item['legend_loc'] = 2
                else:
                    new_item['legend_loc'] = 4

            if not 'legend_title' in item:     # fill legend title
                if item_id in ['h', 'u']:
                    new_item['legend_title'] = dg_label_time
                else:
                    new_item['legend_title'] = ''

            if not 'legend_bbox' in item:      # fill legend bbox_to_anchor
                if item_id in ['h', 'u']:
                    new_item['legend_bbox'] = (1.01, 1.)
                else:
                    new_item['legend_bbox'] = None

            if not 'show_legend' in item:      # fill show legend
                new_item['show_legend'] = True

            # fill data that contains also optional values
            new_item['data'] = []

            ref_num = 1
            for (idx, subitem_data) in enumerate(item_data):
                # copy each item_data and fill with mandatory x and y data

                subitem_data_len = len(subitem_data)

                # resolve label of the plot line
                if (subitem_data_len >=3) and (subitem_data[2] is not None):
                    new_data_subitem_label = subitem_data[2]
                else:
                    if idx == 0:
                        new_data_subitem_label = 'computed'
                    elif idx == 1:
                        new_data_subitem_label = 'measured'
                    else:
                        new_data_subitem_label = 'reference ' + str(ref_num)
                        ref_num += 1

                # resolve plot style of the plot line
                if subitem_data_len >=4:
                    new_data_subitem_style = subitem_data[3]
                else:
                    if item_id in ['h', 'u']:
                        new_data_subitem_style = '-'
                    else:
                        if idx == 0:
                            new_data_subitem_style = '.'
                        else:
                            new_data_subitem_style = 'x'

                # append the newly created data_item
                new_data_item = (subitem_data[0], subitem_data[1],
                                 new_data_subitem_label, new_data_subitem_style)

                new_item['data'].append(new_data_item)

            conformed_plots[item_id].append(new_item)

        return conformed_plots

    def order_plots(narrowed_plots):
        ordered_plots = []
        single_plots  = []

        # first order pairs (and queue singles from pairs)
        for (i1_id, i2_id) in DG_PAIRS:
            i1 = narrowed_plots[i1_id]
            i2 = narrowed_plots[i2_id]

            if i1 and i2:
                ordered_plots.extend(i1)
                ordered_plots.extend(i2)
            elif i1:
                single_plots.extend(i1)
            elif i2:
                single_plots.extend(i2)

            del narrowed_plots[i1_id]
            del narrowed_plots[i2_id]

        # append singles
        ordered_plots.extend(single_plots)
        # append remaning singles that are not part of any pair
        for plots_id in narrowed_plots.values():
            if plots_id: ordered_plots.extend(plots_id)

        return ordered_plots

    def save_text(savedir, plots):
        filename = savedir +'data_as_text.txt'
        fout = open(filename, mode='w', encoding='utf-8')

        for plot in plots:
            plot_id = plot['id']
            ref_num = 1

            for (idx, item) in enumerate(plot['data']):
                if idx == 0:
                    id_suffix = '_comp'
                elif idx == 1:
                    id_suffix = '_meas'
                else:
                    idx_suffix = '_ref' + str(ref_num)
                    ref_num += 1

                fout.write('{:8} = [{}]\n'.format(plot_id + id_suffix + '_x',
                                                  ', '.join(nd2strlist(item[0]))))
                fout.write('{:8} = [{}]\n'.format(plot_id + id_suffix + '_y',
                                                  ', '.join(nd2strlist(item[1]))))
        fout.close()

    def display_plots(ordered_plots):
        nonlocal fignum, save_figures, save_as_text

        if save_figures or save_as_text:
            from shared import get_directories

            if not experiment_info:
                print('Experiment information was not supplied. '
                      'Disabling saving figures.')
                save_figures = save_as_text = False
            else:
                mask = experiment_info['mask']
                if mask:
                    figs_dir_type = 'figs_masks'
                else:
                    figs_dir_type = 'figs_data'

                save_dir = get_directories(figs_dir_type,
                                           experiment_info['exp_id'],
                                           experiment_info['exp_no'],
                                           experiment_info['tube_no'])
                if mask: save_dir += mask + '/'

                if not path.exists(save_dir):
                    makedirs(save_dir)

        if separate_figures:
            images_per_figure = 1
        else:
            images_per_figure = 6

        fignum -= 1
        img_num = 2^20 # high initialization, so that first fig is created

        for plot in ordered_plots:
            plot_id = plot['id']

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
            for plot_data in plot['data']:
                (xdata, ydata, data_label, plot_style) = plot_data
                plt.plot(xdata, ydata, plot_style, label=data_label)

            (xlabel, ylabel) = plot['axes_labels']
            plt.xlabel(xlabel)
            plt.ylabel(ylabel)
            if plot['xscale']:
                plt.xscale(plot['xscale'])
            if plot['yscale']:
                plt.yscale(plot['yscale'])

            if plot['show_legend']:
                plt.legend(borderaxespad=0.0, prop={'family': 'monospace'},
                           loc=plot['legend_loc'], title=plot['legend_title'],
                           bbox_to_anchor=plot['legend_bbox'])

            if save_figures and (img_num == images_per_figure):
                if separate_figures: img_suffix = plot_id
                else: img_suffix = str(fignum)

                plt.savefig(save_dir + 'image-' + img_suffix, dpi=300)

            img_num += 1

        if save_figures and (img_num < images_per_figure):
            plt.savefig(save_dir + 'image-' + str(fignum), dpi=300)

        if save_as_text:
            save_text(save_dir, ordered_plots)

    if plots:
        nplots = narrow_plots(plots)
        oplots = order_plots(nplots)
        display_plots(oplots)

        plt.show(block=False)
        input('Press ENTER to continue...')

        return

    ############################## ORIGINAL CODE ##############################
    if save_figures or save_as_text:
        if not experiment_info:
            print('Experiment information was not supplied. '
                  'Disabling saving figures.')
            save_figures = save_as_text = False
        else:
            mask = experiment_info['mask']
            if mask:
                figs_dir_type = 'figs_masks'
            else:
                figs_dir_type = 'figs_data'

            save_dir = get_directories(figs_dir_type, experiment_info['exp_id'],
                                       experiment_info['exp_no'],
                                       experiment_info['tube_no'])

            if mask: save_dir += mask + '/'

            if not path.exists(save_dir):
                makedirs(save_dir)

    t = [ti/60. for ti in times] # sec -> min

    if has_data(t_ref):
        t_ref = [ti/60. for ti in t_ref] # sec -> min
    else:
        t_ref = t

    if not separate_figures:
        plt.figure(fignum, figsize=(16, 8.5))
        plt.subplots_adjust(wspace=0.15, left=0.06, right=0.85)
    else:
        fignum = 1

    if separate_figures:
        images_per_figure = 1
    else:
        images_per_figure = 6
    img_num = 1

    if has_data(h) or has_data(u):
        if not (has_data(s1) and has_data(s2) and has_data(y)):
            print('draw_graphs error: for ''h'' and/or ''u'' to be displayed '
                  'all ''s1'', ''s2'' and ''y'' have to be set.')
        else:
            x = y2x(y, s1, s2)

            legend_title = "Time [min]"
            legend_data  = ['% 7d' % ti for ti in times]

            if not h is None:
                if separate_figures:
                    plt.figure(fignum)
                else:
                    plt.subplot(321)

                if draw_equilibrium:
                    c_omega = (model.omega ** 2)/model.g/2
                    r0 = model.r0
                    c_eqlib = (r0 + model.l0 + model.fl2) ** 2
                    h_eqlib = c_omega * ((r0 + x) ** 2 - c_eqlib)
                    h_lines = plt.plot(x.transpose(), h.transpose(), '.',
                                       x.transpose(), h_eqlib.transpose(), 'x')

                    h_eq_legend = legend_data + ['h in equilibrium']
                else:
                    h_eq_legend = legend_data
                    h_lines  = plt.plot(x.transpose(), h.transpose(), '.')

                plt.xlabel('Sample length ''L'' [cm]')
                plt.ylabel('Piezometric head ''h'' [cm]')


                if save_figures and separate_figures:
                    plt.savefig(save_dir + 'Image-h', dpi=300)

                add_legend(h_lines, legend_data=h_eq_legend,
                           legend_title=legend_title)

                if separate_figures:
                   add_legend(h_lines, legend_data=h_eq_legend,
                              legend_title=legend_title)
                if save_figures and separate_figures:
                    plt.savefig(save_dir + 'Image-h-leg', dpi=300)

                img_num = 2

            if has_data(u):
                if separate_figures:
                    fignum = fignum + 1
                    plt.figure(fignum)
                else:
                    plt.subplot(3,2,img_num)

                u_lines = plt.plot(x.transpose(), u.transpose(), '.')
                plt.xlabel('Sample length ''L'' [cm]')
                plt.ylabel('Relative saturation ''u''')

                if save_figures and separate_figures:
                    plt.savefig(save_dir + 'Image-u', dpi=300)

                add_legend(u_lines, legend_data=legend_data,
                           legend_title=legend_title)

                if save_figures and separate_figures:
                    plt.savefig(save_dir + 'Image-u-leg', dpi=300)

            img_num = 3

    if (not s1 is None) and all(s1 == 0.0):   s1 = None
    if (not s2 is None) and all(s2 == s2[0]): s2 = None
    if (not mass_in  is None) and all(mass_in  == 0.0): mass_in  = None
    if (not mass_out is None) and all(mass_out == 0.0): mass_out = None

    twins = (((mass_out, mass_out_ref), (mass_in, mass_in_ref)),
             ((GC, GC_ref), (RM, RM_ref)),
             ((s1, s1_ref), (s2, s2_ref)), ((WM, WM_ref), (None, None)))
    ylabels = (('Expelled water [cm]', 'Inflow water [cm]'),
               ('Gravitational center [cm]',
                'Rotational momentum [kg.m.s$^{-1}$]'),
               ('Interface s1 [cm]', 'Interface s2 [cm]'),
               ('Water mass [cm]', None))

    pairs          = []
    pairs_labels   = []
    singles        = []
    singles_labels = []

    # divide all display data into two categories - those to be shown
    # in pairs (like s1 and s2) - aka in the same row; and then the rest
    for (twin, twin_label) in zip(twins, ylabels):
        if twin[0][0] is None or twin[1][0] is None:
            singles.extend(twin)
            singles_labels.extend(twin_label)
        else:
            pairs.extend(twin)
            pairs_labels.extend(twin_label)

    pairs.extend(singles)
    pairs_labels.extend(singles_labels)
    data = pairs
    data_labels = pairs_labels

    nr_measurements     = np.alen(t)
    nr_ref_measurements = np.alen(t_ref)

    for ((ydata, ydata_ref), ydata_label) in zip(data, data_labels):
        if not has_data(ydata): continue

        if img_num > images_per_figure:
            if save_figures:
                plt.savefig(save_dir + ('Image-%i' % fignum), dpi=300)

            img_num = 1
            fignum = fignum + 1

            if separate_figures:
                plt.figure(fignum)
            else:
                plt.figure(fignum, figsize=(16, 8.5))
                plt.subplots_adjust(wspace=0.15, left=0.06, right=0.85)

        if not separate_figures:
            plt.subplot(3,2,img_num)

        # measurements like GC can miss the first (initial) value i.e.
        # values (measurements) are taken only at the end of each centrifuge run
        if np.alen(ydata) == (nr_measurements -1):
            t_plot = t[1:]
        else:
            t_plot = t

        if has_data(ydata_ref):
            if np.alen(ydata_ref) == (nr_ref_measurements -1):
                t_ref_plot = t_ref[1:]
            else:
                t_ref_plot = t_ref

            ydata_lines = plt.plot(t_plot, ydata, '.',
                                   t_ref_plot, ydata_ref, 'x')
            add_legend(ydata_lines, legend_data = ['computed', 'measured'],
                       legend_type='legend', legend_loc=4)
        else:
            plt.plot(t_plot, ydata, '.')

        plt.xlabel('Time [min]')
        plt.ylabel(ydata_label)

        img_num = img_num + 1

    if save_figures:
        plt.savefig(save_dir + ('Image-%i' % fignum), dpi=300)

    plt.show(block=False)

    if save_as_text:
        filename = save_dir +'data_as_text.txt'
        fout = open(filename, mode='w', encoding='utf-8')

        if not h is None:
            fout.write('{:8} = [{}]\n'.format('h', ', '.join(nd2strlist(h))))
        if not u is None:
            fout.write('{:8} = [{}]\n'.format('u', ', '.join(nd2strlist(u))))

        fout.write('{:8} = [{}]\n'.format('t', ', '.join(nd2strlist(t))))
        if not mass_out is None:
            fout.write('{:8} = [{}]\n'.format('mass_out',
                                              ', '.join(nd2strlist(mass_out))))
        if not mass_in is None:
            fout.write('{:8} = [{}]\n'.format('mass_in',
                                              ', '.join(nd2strlist(mass_in))))
        if not GC is None:
            fout.write('{:8} = [{}]\n'.format('GC', ', '.join(nd2strlist(GC))))
        if not RM is None:
            fout.write('{:8} = [{}]\n'.format('RM', ', '.join(nd2strlist(RM))))
        if not WM is None:
            fout.write('{:8} = [{}]\n'.format('WM', ', '.join(nd2strlist(WM))))
        if not s1 is None:
            fout.write('{:8} = [{}]\n'.format('s1', ', '.join(nd2strlist(s1))))
        if not s2 is None:
            fout.write('{:8} = [{}]\n'.format('s2', ', '.join(nd2strlist(s2))))

        if not model.descr == '':
            fout.write('{:8} = "{}"\n'.format('descr', model.descr))
        fout.close()

    input('Press ENTER to continue...')

def disp_inv_results(model, t_inv, inv_params=None,
                     wl_in_inv=None, wl_out_inv=None,
                     gc1_inv=None, rm1_inv=None, cov=None,
                     # results of the direct problem
                     y = None, h_inv = None, u_inv = None,
                     s1_inv = None, s1_ref=None, s2_inv = None, s2_ref=None,
                     wm=None,
                     # other options
                     display_graphs=True, disp_abserror=False,
                     fignum = 1,
                     save_figures=False, separate_figures=False,
                     save_as_text=False, draw_equilibrium=False,
                     show_figures=False, experiment_info=None):

    def print_data(name, data_computed, data_measured):
        name_len = len(name)

        i0 = 0
        in_row = 10
        remaining = np.alen(data_computed)

        if disp_abserror:
            abs_error = np.abs(data_computed - data_measured)
        error = (data_computed - data_measured) / data_measured * 100.

        print('\n')
        while remaining > 0:
            if remaining > in_row:
                disp_items = in_row
            else:
                disp_items = remaining

            print('%s measured: ' % name,
                  disp_items * '% 10.6f' % tuple(data_measured[i0:i0+disp_items]))
            print('%s computed: ' % name,
                  disp_items * '% 10.6f' % tuple(data_computed[i0:i0+disp_items]))
            if disp_abserror:
                print('AbsError: ', name_len * ' ',
                  disp_items * '% 10.6f' % tuple(abs_error[i0:i0+disp_items]))
            print('Error (%):', name_len * ' ',
                  disp_items * '% 10.2f' % tuple(error[i0:i0+disp_items]))

            remaining = remaining - disp_items
            print(108 * '-')
            i0 = i0 + in_row

        print('LSQ error:', np.sum(np.power(data_computed - data_measured, 2)))

    if model.calc_wl_in:
        wl1  = np.asarray(model.get_iterable_value('wl1'), dtype=float)
        print_data('WL_in', wl1_inv, wl1)
    else:
        wl1 = None

    if model.calc_wl_out:
        wl_out_meas  = np.asarray(model.get_iterable_value('wl_out'))
        wl_out_meas  = wl_out_meas.cumsum()
        wl_out_meas[wl_out_meas == 0.0] = 1.0e-10
        print_data('WL_out', wl_out_inv, wl_out_meas)
    else:
        wl_out_meas = None

    if model.calc_gc:
        gc1  = np.asarray(model.get_iterable_value('gc1'), dtype=float)
        print_data('GC', gc1_inv, gc1)
    else:
        gc1 = None
    if model.calc_rm:
        rm1  = np.asarray(model.get_iterable_value('rm1'), dtype=float)
        print_data('RM', rm1_inv, rm1)
    else:
        rm1 = None

    print()
    if inv_params:
        print('Found inverse parameters:')
        for (name, value) in inv_params.items():
            if name == 'ks':
                print('  Ks [cm/s]: {: .8g}'.format(value))
            else:
                print('  {:9}: {: .8g}'.format(name, value))
    if has_data(cov):
        print('Cov:\n', cov)

    if display_graphs:
        from modules.shared.functions import measurements_time

        t_ref = measurements_time(model)

        draw_graphs(t_inv, t_ref=t_ref,
                    mass_in=wl_in_inv, mass_in_ref=wl1,
                    mass_out=wl_out_inv, mass_out_ref=wl_out_meas,
                    GC=gc1_inv, GC_ref=gc1, RM=rm1_inv, RM_ref=rm1,
                    WM=wm, model=model,
                    y=y, h=h_inv, u=u_inv,
                    s1=s1_inv, s1_ref=s1_ref, s2=s2_inv , s2_ref=s2_ref,
                    save_figures=model.save_figures,
                    separate_figures=model.separate_figures,
                    save_as_text=model.save_as_text, draw_equilibrium=False,
                    show_figures=model.show_figures,
                    experiment_info=experiment_info)

def mk_status_item(data_id, data_computed, data_measured = []):
   return {'id': data_id, 'data': (data,computed, data_measured)}

def disp_status(data_plots=None, params=None, cov=None):

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

    if params:
        print('Parameters:')
        for (name, value) in inv_params.items():
            if name == 'ks':
                print('  Ks [cm/s]: {: .8g}'.format(value))
            else:
                print('  {:9}: {: .8g}'.format(name, value))

    if has_data(cov):
        print('Cov:\n', cov)
