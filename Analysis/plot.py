import numpy as np

import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator, FormatStrFormatter

from .low_error import ProcessLowError


def ProcessPlot(cluster, date, app):
    """Controls full process in plotting 2CDs and CMDs.

    CMDs plot V magnitude vs. color, and 2CDs plot R-H vs. color, allowing a
    visual representation of Be candidates.

    Args:
        cluster (str): Cluster from which data is plotted.
        date (str): Date from which data is plotted.
        app (Application): The GUI application object that controls processing.

    """

    file_types = ['', 'lowError']

    for file_type in file_types:
        SinglePlot(cluster, date, app, file_type, '2cd')
        SinglePlot(cluster, date, app, file_type, 'cmd')


def SinglePlot(cluster, date, app, file_type, plot_type):
    """Creates a single plot of given data.

    General configuration of plotting style and other specifications, including
    data, title, labels, and error bars.  It is then output to a file in the
    output directory.

    Args:
        cluster (str): Cluster from which data is plotted.
        date (str): Date from which data is plotted.
        app (Application): The GUI application object that controls processing.
        file_type (str): Distinguishes between regular data or low-error data.
        plot_type (str): Distinguishes between plotting a 2CD or CMD.

    """
    # Setup data set
    path = 'output/' + cluster + '/' + date + '/'

    filename = path + 'phot_scaled_accepted.dat'
    data_in = np.loadtxt(filename, ndmin=2)

    filename = path + 'phot_scaled_rejected.dat'
    data_out = np.loadtxt(filename, ndmin=2)

    filename = path + 'beList_scaled.dat'
    filtered_data = np.loadtxt(filename, ndmin=2)

    if file_type == 'lowError':
        data_in = ProcessLowError(cluster, date, data_in)
        data_out = ProcessLowError(cluster, date, data_out)
        filtered_data = ProcessLowError(cluster, date, filtered_data)

    # Setup plot items
    if plot_type == '2cd':
        y_in = data_in[:, 6] - data_in[:, 8]
        y_err_in = np.sqrt(data_in[:, 7]**2 + data_in[:, 9]**2)
        y_out = data_out[:, 6] - data_out[:, 8]
        y_err_out = np.sqrt(data_out[:, 7]**2 + data_out[:, 9]**2)

        if filtered_data.size > 0:
            be_y = filtered_data[:, 6] - filtered_data[:, 8]
        else:
            be_y = np.array([])

        title = 'R-Halpha vs. B-V'
        y_label = 'R-Halpha'
        output = '2CD_' + file_type + '.png'
    elif plot_type == 'cmd':
        y_in = data_in[:, 4]
        y_err_in = data_in[:, 5]
        y_out = data_out[:, 4]
        y_err_out = data_out[:, 5]

        if filtered_data.size > 0:
            be_y = filtered_data[:, 4]
        else:
            be_y = np.array([])

        title = 'V vs. B-V'
        y_label = 'V'
        output = 'CMD_' + file_type + '.png'

    x_in = data_in[:, 2] - data_in[:, 4]
    x_err_in = np.sqrt(data_in[:, 3]**2 + data_in[:, 5]**2)
    x_out = data_out[:, 2] - data_out[:, 4]
    x_err_out = np.sqrt(data_out[:, 3]**2 + data_out[:, 5]**2)

    if filtered_data.size > 0:
        be_x = filtered_data[:, 2] - filtered_data[:, 4]
    else:
        be_x = np.array([])

    if file_type == 'lowError':
        title += ' (Low Error)'
    x_label = 'B-V'

    # Create plot
    plt.style.use('researchpaper')
    fig, ax = plt.subplots()

    ax.plot(x_out, y_out, 'o', color='#6ba3ff', markersize=11,
            label='Outside cluster')
    ax.plot(x_in, y_in, 'o', color='#3d3d3d', markersize=12,
            label='Inside cluster')

    ax.set_xlabel(x_label)
    ax.set_ylabel(y_label)

    ax.set_xlim([app.B_VMin - 0.1, app.B_VMax + 2.5])
    ax.set_ylim([18.5 - app.A_v, 8.5 - app.A_v])
    ax.set_ylim([18, 7.5])

    ax.xaxis.set_major_locator(MultipleLocator(1))
    ax.xaxis.set_major_formatter(FormatStrFormatter('%d'))
    ax.xaxis.set_minor_locator(MultipleLocator(0.25))

    ax.yaxis.set_major_locator(MultipleLocator(2))
    ax.yaxis.set_major_formatter(FormatStrFormatter('%d'))
    ax.yaxis.set_minor_locator(MultipleLocator(0.5))

    spine_lw = 4
    [ax.spines[axis].set_linewidth(spine_lw)
     for axis in ['top', 'bottom', 'left', 'right']]

    ax.errorbar(x_in, y_in, xerr=x_err_in, yerr=y_err_in, fmt='none',
                ecolor='#8c8c8c', elinewidth=7)
    ax.errorbar(x_out, y_out, xerr=x_err_out, yerr=y_err_out, fmt='none',
                ecolor='#8c8c8c', elinewidth=7)

    # Overplot Be candidates
    ax.plot(be_x, be_y, 'x', color='#ff5151', markersize=15, markeredgewidth=5,
            label='Be Candidates (in and out)')

    # Plot threshold line if 2CD
    if plot_type == '2cd':
        filename = 'standards/' + date + '/' + cluster + \
                   '_aperture_corrections.dat'
        apCorr = np.loadtxt(filename)

        ax.set_ylim([-6.5 + apCorr[2], -4 + apCorr[2]])

        ax.yaxis.set_major_locator(MultipleLocator(1))
        ax.yaxis.set_major_formatter(FormatStrFormatter('%d'))
        ax.yaxis.set_minor_locator(MultipleLocator(0.25))

        filename = 'output/' + cluster + '/' + date + \
                   '/thresholds' + '.dat'
        thresholds = np.loadtxt(filename)
        if app.threshold_type == 'Constant':
            file = thresholds[0]
        elif app.threshold_type == 'Linear':
            file = thresholds[1]
        slope = file[0]
        intercept = file[1]

        # linex = np.array([app.B_VMin, app.B_VMax])
        linex = np.array([app.B_VMin - 0.1, app.B_VMax + 2.5])
        liney = slope * linex + intercept
        ax.plot(linex, liney, '--', color='#ff5151', label='Be Threshold',
                linewidth=6)

    ax.legend()

    # Output
    filename = 'output/' + cluster + '/' + date + '/plots/' + output
    fig.savefig(filename)

    plt.close('all')
