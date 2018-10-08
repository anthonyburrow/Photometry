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

        ColorMagnitudeDiagram(cluster, date, app, data, filtered_data, file_type)
        TwoColorDiagram(cluster, date, app, data, filtered_data, file_type)


def ColorMagnitudeDiagram(cluster, date, app, data, filtered_data, file_type):
    """Specifies the plotting of a color magnitude diagram.

    Calls to plot a color magnitude diagram with V vs. B-V axes.

    """
    print("  Creating color-magnitude diagram...")

    # Plot full data
    B_V = data[:, 2] - data[:, 4]
    B_Verr = np.sqrt(data[:, 3]**2 + data[:, 5]**2)

    try:
        filtered_B_V = filtered_data[:, 2] - filtered_data[:, 4]
        filtered_V = filtered_data[:, 4]
    except IndexError:
        filtered_B_V = []
        filtered_V = []

    title = 'V vs. B-V'
    if file_type == '_lowError':
        title += ' (Low Error)'
    SinglePlot(cluster, date, app, B_V, data[:, 4], title, 'B-V', 'V',
               B_Verr, data[:, 5], filtered_B_V, filtered_V,
               'CMD_' + app.phot_type + file_type + '.png')


def TwoColorDiagram(cluster, date, app, data, filtered_data, file_type):
    """Specifies the plotting of a two-color diagram.

    Calls to plot a color magnitude diagram with R-H vs. B-V axes.

    """
    print("  Creating two-color diagram...")

    # Plot full data
    B_V = data[:, 2] - data[:, 4]
    B_Verr = np.sqrt(data[:, 3]**2 + data[:, 5]**2)
    R_H = data[:, 6] - data[:, 8]
    R_Herr = np.sqrt(data[:, 7]**2 + data[:, 9]**2)

    try:
        filtered_B_V = filtered_data[:, 2] - filtered_data[:, 4]
        filtered_R_H = filtered_data[:, 6] - filtered_data[:, 8]
    except IndexError:
        filtered_B_V = []
        filtered_R_H = []

    title = 'R-Halpha vs. B-V'
    if file_type == '_lowError':
        title += ' (Low Error)'
    SinglePlot(cluster, date, app, B_V, R_H, title, 'B-V', 'R-Halpha',
               B_Verr, R_Herr, filtered_B_V, filtered_R_H,
               '2CD_' + app.phot_type + file_type)


def SinglePlot(cluster, date, app, x, y, title, x_label, y_label, x_err, y_err, be_x, be_y, output):
    """Creates a single plot of given data.

    General configuration of plotting style and other specifications, including data,
    title, labels, and error bars.  It is then output to a file in the output
    directory.

    Args:
            x: The x-axis array to be plotted.
            y: The y-axis array to be plotted.
            title: The title of the plot.
            x_label: The label of the x-axis.
            y_label: The label of the y-axis.
            x_err: The x-axis error information.
            y_err: The y-axis error information.
            be_x: The x-axis Be-filtered data.
            be_y:The y-axis Be-filtered data.
            output: The filename desired for the plot.

    """
    # plt.style.use('ggplot')

    plt.figure(figsize=(12, 9))

    # Plot main data
    plt.plot(x, y, 'o', color='#3f3f3f', markersize=12)
    # plt.title(title)
    plt.xlabel(x_label, fontsize=36)
    plt.ylabel(y_label, fontsize=36)

    plt.axes().xaxis.set_major_locator(MultipleLocator(1))
    plt.axes().xaxis.set_major_formatter(FormatStrFormatter('%d'))
    plt.axes().xaxis.set_minor_locator(MultipleLocator(0.25))

    plt.axes().yaxis.set_major_locator(MultipleLocator(2))
    plt.axes().yaxis.set_major_formatter(FormatStrFormatter('%d'))
    plt.axes().yaxis.set_minor_locator(MultipleLocator(0.5))

    plt.axes().tick_params('both', length=12, width=4, which='major', top=True, right=True, direction='in', pad=6, labelsize=30)
    plt.axes().tick_params('both', length=8, width=3, which='minor', top=True, right=True, direction='in')

    plt.axes().spines['top'].set_linewidth(4)
    plt.axes().spines['right'].set_linewidth(4)
    plt.axes().spines['bottom'].set_linewidth(4)
    plt.axes().spines['left'].set_linewidth(4)

    plt.xlim([app.B_VMin - 0.1, app.B_VMax + 2.5])
    plt.ylim([18.5 - app.A_v, 8.5 - app.A_v])
    plt.errorbar(x, y, xerr=x_err, yerr=y_err, fmt='none', ecolor='#50a0e5', elinewidth=7)

    # Overplot Be candidates
    plt.plot(be_x, be_y, 'x', color='#ff5151', markersize=15, markeredgewidth=5, label='Be Candidates')

    # Plot threshold line if 2CD
    if output == '2CD_' + app.phot_type or output == '2CD_' + app.phot_type + '_lowError':
        # plt.ylim([-5 - self.app.A_r, -1 - self.app.A_r])
        try:
            apCorr = np.loadtxt('standards/' + date + '/' + cluster + '_aperture_corrections.dat')
        except IOError:
            pass

        plt.ylim([-6.5 + apCorr[2], -4 + apCorr[2]])

        plt.axes().yaxis.set_major_locator(MultipleLocator(1))
        plt.axes().yaxis.set_minor_locator(MultipleLocator(0.25))

        filename = 'output/' + cluster + '/' + date + '/thresholds_' + app.phot_type + '.dat'
        try:
            thresholds = np.loadtxt(filename)
            if app.threshold_type == 'Constant':
                file = thresholds[0]
            elif app.threshold_type == 'Linear':
                file = thresholds[1]
            slope = file[0]
            intercept = file[1]

            linex = np.array([app.B_VMin, app.B_VMax])
            liney = slope * linex + intercept
            plt.plot(linex, liney, '--', color='#ff5151', label='Be Threshold', linewidth=6)
        except IOError:
            print("\nNote: Thresholds have not been calculated or written to file yet and will not be displayed.")

    plt.legend(fontsize=28)
    plt.tight_layout()

    # Output
    filename = 'output/' + cluster + '/' + date + '/plots/' + output
    plt.savefig(filename)

    plt.clf()
