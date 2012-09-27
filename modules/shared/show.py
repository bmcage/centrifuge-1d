import matplotlib.pyplot as plt
import numpy as np
import pickle
from const import FIGS_DIR
from os import makedirs, path
from shared import get_directories


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

def add_dplotline(dplot, xdata, ydata, label=None, line_opts='.'):
    if not (has_data(xdata) and has_data(ydata)): return None

    if label is None:
        label = 'reference ' + str(dplot['_ref_num'])
        dplot['_ref_num'] += 1

    item = (xdata, ydata, label, line_opts)

    dplot['data'].append(item)

    return item

def make_dplot(dplot_id, legend_loc=None, show_legend=None, legend_bbox=None,
            legend_title=dg_label_time, xscale=None, yscale=None,
            axes_labels=None):

    if axes_labels is None: # try to set default value
        if dplot_id in DG_AXES_LABELS:
            axes_labels = DG_AXES_LABELS[dplot_id]

    if ((dplot_id in ['h', 'u']) and (legend_bbox is None)
        and (legend_loc is None)):
        legend_bbox = (1.01, 1.)
        legend_loc = 2

    if legend_loc is None:
        legend_loc = 4

    dplot = {'id': dplot_id, 'data': [], 'axes_labels': axes_labels,
             'legend_title': legend_title, 'legend_bbox': legend_bbox,
             'show_legend': show_legend, 'xscale': xscale, 'yscale': yscale,
             'legend_loc': legend_loc, '_ref_num': 1}

    return dplot

def display_dplots(dplots, save_figures=False, separate_figures=False,
                   save_as_text=False,show_figures=False, experiment_info=None,
                   fignum = 1):

    if not dplots: return

    def order_dplots(dplots):
        ordered_dplots = []
        single_dplots  = []

        dplots_bucket = {}
        for dplot in dplots:
            dplot_id = dplot['id']

            if not dplot['data']: continue

            if dplot_id in dplots_bucket:
                dplots_bucket[dplot_id].append(dplot)
            else:
                dplots_bucket[dplot_id] = [dplot]

        # first order pairs (and queue singles from pairs)
        for (i1_id, i2_id) in DG_PAIRS:
            i1p = (i1_id in dplots_bucket)
            i2p = (i2_id in dplots_bucket)

            if i1p and i2p:
                ordered_dplots.extend(dplots_bucket[i1_id])
                ordered_dplots.extend(dplots_bucket[i2_id])
            elif i1p:
                single_dplots.extend(dplots_bucket[i1_id])
            elif i2p:
                single_dplots.extend(dplots_bucket[i2_id])

            if i1p: del dplots_bucket[i1_id]
            if i2p: del dplots_bucket[i2_id]

        # append singles
        ordered_dplots.extend(single_dplots)
        # append remaning singles that are not part of any pair
        for remaining_dplots in dplots_bucket.values():
            ordered_dplots.extend(remaining_dplots)

        return ordered_dplots

    def save_text(savedir, dplots):
        filename = savedir +'data_as_text.txt'
        fout = open(filename, mode='w', encoding='utf-8')

        for plot in dplots:
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

    def show_dplots(ordered_dplots):
        nonlocal fignum, save_figures, save_as_text

        if save_figures or save_as_text:
            from shared import get_directories

            if not experiment_info:
                print('Experiment information was not supplied. '
                      'Disabling saving figures.')
                save_figures = save_as_text = False
            else:
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
                (xdata, ydata, data_label, plot_style) = dplot_data
                plt.plot(xdata, ydata, plot_style)
                if type(data_label) == str:
                    plot_labels.append(data_label)
                else:
                    plot_labels.extend(data_label)

            (xlabel, ylabel) = dplot['axes_labels']
            plt.xlabel(xlabel)
            plt.ylabel(ylabel)
            if dplot['xscale']:
                plt.xscale(dplot['xscale'])
            if dplot['yscale']:
                plt.yscale(dplot['yscale'])

            show_legend = dplot['show_legend']
            if show_legend is None:
                if (len(dplot['data']) > 1) or (np.ndim(ydata) > 1):
                 show_legend = True
            if show_legend:
                plt.legend(plot_labels, borderaxespad=0.0, prop={'family': 'monospace'},
                           loc=dplot['legend_loc'], title=dplot['legend_title'],
                           bbox_to_anchor=dplot['legend_bbox'])

            if save_figures and (img_num == images_per_figure):
                if separate_figures: img_suffix = dplot_id
                else: img_suffix = str(fignum)

                plt.savefig(save_dir + 'image-' + img_suffix, dpi=300)

            img_num += 1

        if save_figures and (img_num < images_per_figure):
            plt.savefig(save_dir + 'image-' + str(fignum), dpi=300)

        if save_as_text:
            save_text(save_dir, ordered_dplots)

    if type(dplots) == dict: # we have single dplot item
        dplots = (dplots, )
    ord_dplots = order_dplots(dplots)
    show_dplots(ord_dplots)

    plt.show(block=False)
    input('Press ENTER to continue...')

def draw_graphs(times, t_ref = None, y = None, h = None, u = None,
                s1 = None, s1_ref=None, s2 = None, s2_ref=None,
                mass_out = None, mass_out_ref = None,
                mass_in = None, mass_in_ref = None,
                GC = None, GC_ref = None,
                RM = None,  RM_ref = None, WM = None, WM_ref = None,
                fignum = 1, save_figures=False, separate_figures=False,
                save_as_text=False, draw_equilibrium=False,
                show_figures=False, experiment_info=None,
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

    if save_figures or save_as_text:
        if not experiment_info:
            print('Experiment information was not supplied. '
                  'Disabling saving figures.')
            save_figures = save_as_text = False
        else:
            mask = experiment_info['mask']
            if experiment_info['mask']:
                figs_dir_type = 'mask'
            else:
                figs_dir_type = 'data'

            save_dir = get_directories('figs', figs_dir_type,
                                       experiment_info)

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
   return {'id': data_id, 'data': (data_computed, data_measured)}

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
                print('M', value_measured, '\nC', value_computed)
                abs_error = np.abs(value_computed - value_measured)

            compare_data(plot_id, value_computed, value_measured,
                         rel_error, abs_error)

    if params:
        print('Parameters:')
        for (name, value) in params.items():
            if name == 'ks':
                print('  Ks [cm/s]: {: .8g}'.format(value))
            else:
                print('  {:9}: {: .8g}'.format(name, value))

    if has_data(cov):
        print('Cov:\n', cov)

# ResultData: hold the data of the computation
# Internal structure: dictionary of types: {line_type: line_type_data} where
# line_type is (by default) one of ['computed', 'measured', 'references'] and
# line_type_data is a dictionary of types: {line_id: line_data}, with
# line_id the ID of the line and line_data a dict of types
# {data_type: (xdata, ydata)}, where data_type is (by default) one of
# ['h', 'u', 'GC', 'RM', 'WM', 's1','s2'], xdata is the x-axis coordinate and
# ydata is the y-axis coordinate (computed value)
class ResultsData():
    def __init__(self, extract_data_fn, model, referencing_parameters=[],
                 measurements=None):
        # add computed data
        (flag, value) = extract_data_fn(model)
        if not flag:
            print('Computation was not successfull. No data will be saved.')
            self._data = None
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

        self._data = data
    def has_data(self):
        return not self._data is None

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
            exit(1)

        with open(filename, 'rb') as f:
            pickle.load(self._data, f)
    def get_datatypes(self):
         # line_type 'computed' with ID 'computed' has to be present
        return self._data['computed']['computed'].keys()
    def get_linetypes(self):
        return self._data.keys()

    def get_linedata(self, line_type, line_id):
        data = self._data
        if (line_type in data) and (line_id in data[line_type]):
            return data[line_type][line_id]
        else:
            return None

    def iterate(self):
        for (data_type, lines) in self._data.items():
            for (line_id, line_data) in lines.items():
                yield (data_type, line_id, line_data)

PLOTSTYLE_ININAME = 'plotstyle.ini'

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

            plotstyle_files = \
              filter_existing(prefix_with_paths(PLOTSTYLE_ININAME, search_dirs))

            return plotstyle_files

        def read_plotstyle_cfg(filename):
            result = {}

            # read config files
            parser   = configparser.ConfigParser()
            try:
                read_files = parser.read(fname)
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
        plotstyle_filenames = get_filenames(experiment_info)

        if not plotstyle_filenames: return {}

        plot_cfg = read_plotcfg(plotstyle_filenames[0])

        for fname in plotstyle_filenames[1:]:
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
                                              'legend_bbox', 'legend_loc']}
        user_styles = self._get_value('datasets')
        if user_styles: dplot_styles.update(user_styles)

         # try to set default value
        if (dplot_styles['axes_labels'] is None) and (dplot_id in DG_AXES_LABELS):
           dplot_styles['axes_labels'] = DG_AXES_LABELS[dplot_id]

        if ((dplot_id in ['h', 'u']) and (dplot_styles['legend_bbox'] is None)
            and (dplot_styles['legend_loc'] is None)):
            dplot_styles['legend_bbox'] = (1.01, 1.)
            dplot_styles['legend_loc'] = 2

        if (dplot_id in ['h', 'u']) and (dplot_styles['legend_title'] is None):
            dplot_styles['legend_title'] = 'Time [min]'

        if ((dplot_id == 'h') and (dplot_styles['show_legend'] is None)
            and (not self.get_display_options()['separate_figures'])):
            dplot_styles['show_legend'] = False

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

        user_styles = self._get_value('lines')
        if user_styles and (line_id in user_styles):
            line_styles.update(user_styles[line_id])

        return line_styles

    def get_linestyle(self, line_type, line_id, data_types):
        return self._mk_linestyle(line_type, line_id, data_types)

class DPlots():
    def __init__(self, data, experiment_info):
        self._data = data
        self._display = data.has_data()
        if not self._display:
            print('No data is provided. Nothing to display.')
            return

        self._experiment_info = experiment_info

        plot_styles = PlotStyles(experiment_info)
        self._plot_styles = plot_styles

        data_types = data.get_datatypes()
        self._dplots_bucket = self._make(data_types, plot_styles)

        for (line_type, line_id, line_data) in data.iterate():
            self._add_plotline(line_type, line_id, line_data, data_types)

    def _make(self, data_types, plot_styles):
        dplots_bucket = {dtype: {'id': dtype, 'data': [],
                                  'styles': plot_styles.get_dplotstyles(dtype)}
                         for dtype in data_types}

        return dplots_bucket

    def _add_plotline(self, line_type, line_id, line_data, data_types):
        bucket = self._dplots_bucket

        line_styles = self._plot_styles.get_linestyle(line_type, line_id,
                                                      data_types)
        if 'label' in line_styles:
            label = line_styles['label']
        else:
            label = line_id

        for (data_type, data_value) in line_data.items():
            if not data_value: continue

            (xdata, ydata) = (data_value[0], data_value[1])

            if not (has_data(xdata) and has_data(ydata)): continue

            if ((data_type in ['h', 'u']) and (len(data_value) > 2)
                and has_data(data_value[2]) and (not 'label' in line_styles)):
                ilabel = ['% 6d' % (ti/60.) for ti in data_value[2]]
            else:
                ilabel = label

            item = (xdata, ydata, ilabel, line_styles[data_type])

            bucket[data_type]['data'].append(item)

    def show_status(self):
        data = self._data
        status_items = []

        measurements = data.get_linedata('measured', 'measured')
        computed     = data.get_linedata('computed', 'computed')
        print('HHH', measurements, computed)
        for (key, m_value) in measurements.items():
            if m_value[1] is None: continue

            if key in computed:
                c_value = computed[key][1]

            if c_value is not None:
                status_items.append(mk_status_item(key, c_value[:-1],
                                                   m_value[1]))

        if status_items:
            disp_status(status_items)

    def display(self, fignum = 1):
        if not self._display: return

        def order_dplots(dplots_bucket):
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

        def show_dplots(ordered_dplots):
            nonlocal fignum

            disp_opts = self._plot_styles.get_display_options()

            separate_figures = disp_opts['separate_figures']
            save_figures     = disp_opts['save_figures']

            if save_figures:
                experiment_info = self._experiment_info
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
                    (xdata, ydata, data_label, plot_style) = dplot_data
                    plt.plot(xdata, ydata, plot_style)
                    if type(data_label) == str:
                        plot_labels.append(data_label)
                    else:
                        plot_labels.extend(data_label)

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

        # function body
        ord_dplots = order_dplots(self._dplots_bucket)
        show_dplots(ord_dplots)

        plt.show(block=False)
        input('Press ENTER to continue...')
