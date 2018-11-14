import numpy as np

import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator, FormatStrFormatter


def ProcessPlot(cluster, date, app):
    """Reads finalized photometry data to be plotted.

    Converts read data file to a multidimensional array used to plot data.
    """
    path = 'output/' + cluster + '/' + date + '/'

    file_types = ['', '_lowError']

    for file_type in file_types:
        filename = path + 'phot_scaled_' + app.phot_type + file_type + '.dat'
        data = np.loadtxt(filename, ndmin=2)

        filename = path + 'beList_scaled_' + app.phot_type + file_type + '.dat'
        filtered_data = np.loadtxt(filename, ndmin=2)

        SinglePlot(cluster, date, app, data, filtered_data, file_type, '2cd')
        SinglePlot(cluster, date, app, data, filtered_data, file_type, 'cmd')


def SinglePlot(cluster, date, app, data, filtered_data, file_type, plot_type):
    """Creates a single plot of given data.

    General configuration of plotting style and other specifications, including
    data, title, labels, and error bars.  It is then output to a file in the
    output directory.

    Args:

    """
    # Setup plot items
    if plot_type == '2cd':
        y = data[:, 6] - data[:, 8]
        y_err = np.sqrt(data[:, 7]**2 + data[:, 9]**2)
        be_y = filtered_data[:, 6] - filtered_data[:, 8]

        title = 'R-Halpha vs. B-V'
        y_label = 'R-Halpha'
        output = '2CD_' + app.phot_type + file_type + '.png'
    elif plot_type == 'cmd':
        y = data[:, 4]
        y_err = data[:, 5]
        be_y = filtered_data[:, 4]

        title = 'V vs. B-V'
        y_label = 'V'
        output = 'CMD_' + app.phot_type + file_type + '.png'

    # try:
    #     filtered_B_V = filtered_data[:, 2] - filtered_data[:, 4]
    #     filtered_V = filtered_data[:, 4]
    # except IndexError:
    #     filtered_B_V = []
    #     filtered_V = []

    x = data[:, 2] - data[:, 4]
    x_err = np.sqrt(data[:, 3]**2 + data[:, 5]**2)
    be_x = filtered_data[:, 2] - filtered_data[:, 4]

    if file_type == '_lowError':
        title += ' (Low Error)'
    x_label = 'B-V'

    # Create plot
    plt.style.use('researchpaper')
    fig, ax = plt.subplots()

    # Plot main data
    ax.plot(x, y, 'o', color='#3f3f3f', markersize=12)
    # plt.set_title(title)
    ax.set_xlabel(x_label)
    ax.set_ylabel(y_label)

    ax.set_xlim([app.B_VMin - 0.1, app.B_VMax + 2.5])
    ax.set_ylim([18.5 - app.A_v, 8.5 - app.A_v])

    ax.xaxis.set_major_locator(MultipleLocator(1))
    ax.xaxis.set_major_formatter(FormatStrFormatter('%d'))
    ax.xaxis.set_minor_locator(MultipleLocator(0.25))

    ax.yaxis.set_major_locator(MultipleLocator(2))
    ax.yaxis.set_major_formatter(FormatStrFormatter('%d'))
    ax.yaxis.set_minor_locator(MultipleLocator(0.5))

    spine_lw = 4
    [ax.spines[axis].set_linewidth(spine_lw)
     for axis in ['top', 'bottom', 'left', 'right']]

    ax.errorbar(x, y, xerr=x_err, yerr=y_err, fmt='none', ecolor='#50a0e5',
                elinewidth=7)

    # Overplot Be candidates
    ax.plot(be_x, be_y, 'x', color='#ff5151', markersize=15, markeredgewidth=5,
            label='Be Candidates')

    # Plot threshold line if 2CD
    if plot_type == '2cd':
        apCorr = np.loadtxt('standards/' + date + '/' + cluster +
                            '_aperture_corrections.dat')

        ax.set_ylim([-6.5 + apCorr[2], -4 + apCorr[2]])

        ax.yaxis.set_major_locator(MultipleLocator(1))
        ax.yaxis.set_major_formatter(FormatStrFormatter('%d'))
        ax.yaxis.set_minor_locator(MultipleLocator(0.25))

        filename = 'output/' + cluster + '/' + date + \
                   '/thresholds_' + app.phot_type + '.dat'
        thresholds = np.loadtxt(filename)
        if app.threshold_type == 'Constant':
            file = thresholds[0]
        elif app.threshold_type == 'Linear':
            file = thresholds[1]
        slope = file[0]
        intercept = file[1]

        linex = np.array([app.B_VMin, app.B_VMax])
        liney = slope * linex + intercept
        ax.plot(linex, liney, '--', color='#ff5151', label='Be Threshold',
                linewidth=6)

    ax.legend()

    # Output
    filename = 'output/' + cluster + '/' + date + '/plots/' + output
    fig.savefig(filename)

    # plt.cla()
    # plt.clf()
