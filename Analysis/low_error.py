import numpy as np

import os.path

import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator, FormatStrFormatter


def ProcessLowError(path, file, app):
    # Get max error for observation date
    filename = path + 'phot_' + app.phot_type + '.dat'
    data = np.loadtxt(filename)
    max_error = CalcMaxError(data)

    # Filter given data based on max error
    filename = path + file + '_' + app.phot_type + '.dat'
    data = np.loadtxt(filename)
    lowError_data = GetLowErrorData(data, max_error)

    # Write to file
    filename = filename[:-4] + '_lowError.dat'
    with open(filename, 'w') as F:
        np.savetxt(F, lowError_data, fmt='%.3f')

    PlotLowError(data, max_error, path, app)


def CalcMaxError(data):
    R_Herr = np.sqrt(data[:, 7]**2 + data[:, 9]**2)

    # Calculate max error by std. of error
    print("\nCalculating the max error for low-error data...")

    std = np.sqrt(1 / (len(R_Herr) - 1)) * np.sqrt(np.sum(R_Herr**2))
    max_error = np.sqrt(2) * std

    print("  Max error found to be: ", "%.3f" % max_error)

    return max_error


def GetLowErrorData(data, max_error):
    lowError_data = []

    R_Herr = np.sqrt(data[:, 7]**2 + data[:, 9]**2)
    for i in range(0, len(data)):
        if R_Herr[i] < max_error:
            lowError_data.append(data[i])

    return lowError_data


def PlotLowError(data, max_error, path, app):
    x = data[:, 6] - data[:, 8]
    y = np.sqrt(data[:, 7]**2 + data[:, 9]**2)
    std = max_error / np.sqrt(2)

    # plt.style.use('ggplot')

    plt.figure(figsize=(12, 9))

    # Plot main data
    plt.plot(x, y, 'o', color='#3f3f3f', markersize=12)
    # plt.title("R-H Err vs. R-H")
    plt.xlabel("R-H", fontsize=36)
    plt.ylabel("R-H Err", fontsize=36)
    plt.xlim([max(x), min(x)])
    plt.ylim([min(y) - 0.01, max(y) + 0.01])

    plt.axes().xaxis.set_major_locator(MultipleLocator(2))
    plt.axes().xaxis.set_major_formatter(FormatStrFormatter('%d'))
    plt.axes().xaxis.set_minor_locator(MultipleLocator(0.5))

    plt.axes().yaxis.set_major_locator(MultipleLocator(0.02))
    plt.axes().yaxis.set_major_formatter(FormatStrFormatter('%.2f'))
    plt.axes().yaxis.set_minor_locator(MultipleLocator(0.005))

    plt.axes().tick_params('both', length=12, width=4, which='major', top=True, right=True, direction='in', pad=6, labelsize=30)
    plt.axes().tick_params('both', length=8, width=3, which='minor', top=True, right=True, direction='in')

    plt.axes().spines['top'].set_linewidth(4)
    plt.axes().spines['right'].set_linewidth(4)
    plt.axes().spines['bottom'].set_linewidth(4)
    plt.axes().spines['left'].set_linewidth(4)

    plt.hlines(std, min(x), max(x), linestyles='dashed', label='Standard Error')

    plt.legend(fontsize=28)
    plt.tight_layout()

    # Output
    filename = path + 'plots/'
    if not os.path.exists(filename):
        os.makedirs(filename)

    filename = path + 'plots/magErr_vs_mag_' + app.phot_type + '.png'
    plt.savefig(filename)

    plt.clf()
