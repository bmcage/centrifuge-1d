import matplotlib.pyplot as plt
from numpy import alen, empty

def y2x(y, s1, s2):
    s1_len = alen(s1)
    if s1_len != alen(s2):
        print('Interfaces array ''s1'' and ''s2'' have to be of the same'
              'lenght. Cannot proceed.')
        exit(1)
    x = empty([s1_len, len(y)], float)

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

SAVE_PATH = '/home/archetyp/projekty/lyx/geosound/'

def draw_graphs(t, y = None, s1 = None, s2 = None, h = None, u = None,
                mass_out = None, mass_in = None,
                GC = None, RM = None,  WM = None,
                fignum = 1, save_figures=False, save_figures_separately=False):

    twins = ((mass_out, mass_in), (GC, RM), (s1, s2), (WM, None))
    ylabels = (('Exspelled water [cm]', 'Inflow water [cm]'),
               ('Gravitational center [cm]',
                'Rotational momentum [kg.m.s$^{-1}$]'),
               ('Interface s1 [cm]', 'Interface s2 [cm]'),
               ('Water mass [cm]', None))

    plt.figure(fignum, figsize=(16, 8.5))

    plt.subplots_adjust(wspace=0.15, left=0.06, right=0.85)

    row = -1

    if (not h is None) and (not u is None):
        if (s1 is None) or (s2 is None) or (y is None):
            print('draw_graphs error: for ''h'' and/or ''u'' to be displayed '
                  'all ''s1'', ''s2'' and ''y'' have to be set.')
        else:
            legend_data = []
            for i in range(len(t)):
                legend_data.append('t =%7d' % t[i])

            x = y2x(y, s1, s2)

            if not h is None:
                plt.subplot(321)
                plt.plot(x.transpose(), h.transpose(), '.')
                plt.xlabel('Rotational axis distance ''r'' [cm]')
                plt.ylabel('Piezometric head ''h'' [cm]')

            if not u is None:
                plt.subplot(322)
                plt.plot(x.transpose(), u.transpose(), '.')
                plt.xlabel('Rotational axis distance ''r'' [cm]')
                plt.ylabel('Relative saturation ''u''')
                plt.legend(legend_data, bbox_to_anchor=(1.02, 1.), loc=2,
                           borderaxespad=0.0)

                row = 0

    for ((v1, v2), (ylabel1, ylabel2)) in zip(twins, ylabels):
        if (not v1 is None) or (not v2 is None):
            if row == 2:
                if save_figures:
                    fname = SAVE_PATH + ('Image-%i' % fignum)
                    plt.savefig(fname, dpi=300)
                fignum = fignum + 1
                row = 0

                plt.figure(fignum, figsize=(16, 8.5))
                plt.subplots_adjust(wspace=0.15, left=0.06, right=0.85)
            else:
                row = row + 1
        else:
            continue

        column = 1

        if not v1 is None:
            plt.subplot(3,2,2*row + column)
            plt.plot(t, v1, '.')
            plt.xlabel('Time [s]')
            plt.ylabel(ylabel1)

            column = 2

        if not v2 is None:
            plt.subplot(3,2,2*row + column)
            plt.plot(t, v2, '.')
            plt.xlabel('Time [s]')
            plt.ylabel(ylabel2)

    plt.show(block=False)
    if save_figures:
        fname = SAVE_PATH + ('Image-%i' % fignum)
        plt.savefig(fname, dpi=300)

    input('Press ENTER to continue...')
