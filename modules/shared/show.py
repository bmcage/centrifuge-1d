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

def draw_graphs(t, y = None, s1 = None, s2 = None, h = None, u = None,
                mass_out = None, mass_in = None,
                GC = None, RM = None,  WM = None,
                fignum = 1):

    twins = ((mass_out, mass_in), (GC, RM), (s1, s2), (WM, None))
    ylabels = (('Outspelled water [$cm$]', 'Inflow water [$cm$]'),
               ('Gravitational center [$cm$]',
                'Rotational momentum [$kg.m.s^{-1}$]'),
               ('Interface s1 [$cm$]', 'Interface s2 [$cm$]'),
               ('Water mass [$cm$]', None))

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
                plt.xlabel('Rotational axis distance ''r'' [$cm$]')
                plt.ylabel('Piezometric head ''h'' [$cm$]')

            if not u is None:
                plt.subplot(322)
                plt.plot(x.transpose(), u.transpose(), '.')
                plt.xlabel('Rotational axis distance ''r'' [$cm$]')
                plt.ylabel('Relative saturation ''u''')
                plt.legend(legend_data, bbox_to_anchor=(1.02, 1.), loc=2,
                           borderaxespad=0.0)

                row = 0

    for ((v1, v2), (ylabel1, ylabel2)) in zip(twins, ylabels):
        if (not v1 is None) or (not v2 is None):
            if row == 2:
                print('in')
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
            plt.xlabel('Time [$s$]')
            plt.ylabel(ylabel1)

            column = 2

        if not v2 is None:
            plt.subplot(3,2,2*row + column)
            plt.plot(t, v2, '.')
            plt.xlabel('Time [$s$]')
            plt.ylabel(ylabel2)

    plt.show()
