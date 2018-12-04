import numpy as np

import os.path


def ProcessLowError(cluster, date, data):
    """Controls full process for low error calculations.

    Low error data is calculated using the standard deviation of the errors
    of a given data set.  Those targets with higher error are then filtered
    out of the set.

    Args:
        data (numpy.array): Data to filter for low error.

    """
    # Get max error for observation date
    filename = 'output/' + cluster + '/' + date + '/phot.dat'
    phot = np.loadtxt(filename, ndmin=2)
    max_error = CalcMaxError(phot)

    # Filter given data based on max error
    if data.size > 0:
        lowError_data = GetLowErrorData(data, max_error)
    else:
        lowError_data = np.array([])

    # PlotLowError(data, max_error, path)

    return lowError_data


def CalcMaxError(data):
    """Calculates the error threshold that defines maximum error.

    Args:
        data (list): Data set for which max error is calculated.

    Returns:
        float: The maximum error threshold.

    """
    R_Herr = np.sqrt(data[:, 7]**2 + data[:, 9]**2)

    # Calculate max error by std. of error
    std = np.sqrt(1 / (len(R_Herr) - 1)) * np.sqrt(np.sum(R_Herr**2))
    max_error = np.sqrt(2) * std

    return max_error


def GetLowErrorData(data, max_error):
    """Filters data set for high error targets.

    Args:
        data (list): Data set to be filtered.
        max_error (float): The calculated max error threshold.

    Returns:
        list: Newly filtered data set.

    """
    lowError_data = []

    for target in data:
        R_Herr = np.sqrt(target[7]**2 + target[9]**2)
        if R_Herr < max_error:
            lowError_data.append(target)

    return np.array(lowError_data)


def PlotLowError(data, max_error, path):
    import matplotlib.pyplot as plt
    from matplotlib.ticker import MultipleLocator, FormatStrFormatter

    """Plots the error distribution.

    Args:
        data (list): Data set being processed.
        max_error (float): The calculated max error threshold.
        path (str): The directory/path that holds output photometry files.

    """
    x = data[:, 6] - data[:, 8]
    y = np.sqrt(data[:, 7]**2 + data[:, 9]**2)
    std = max_error / np.sqrt(2)

    plt.style.use('researchpaper')
    fig, ax = plt.subplots()

    # Plot main data
    ax.plot(x, y, 'o', color='#3f3f3f', markersize=12)

    ax.set_xlabel("R-H")
    ax.set_ylabel("R-H Err")
    ax.set_xlim([max(x), min(x)])
    ax.set_ylim([min(y) - 0.01, max(y) + 0.01])

    ax.xaxis.set_major_locator(MultipleLocator(2))
    ax.xaxis.set_major_formatter(FormatStrFormatter('%d'))
    ax.xaxis.set_minor_locator(MultipleLocator(0.5))

    ax.yaxis.set_major_locator(MultipleLocator(0.02))
    ax.yaxis.set_major_formatter(FormatStrFormatter('%.2f'))
    ax.yaxis.set_minor_locator(MultipleLocator(0.005))

    spine_lw = 4
    [ax.spines[axis].set_linewidth(spine_lw)
     for axis in ['top', 'bottom', 'left', 'right']]

    ax.hlines(std, min(x), max(x), linestyles='dashed', label='Standard Error')

    # Output
    filename = path + 'plots/'
    if not os.path.exists(filename):
        os.makedirs(filename)

    filename = path + 'plots/magErr_vs_mag.png'
    fig.savefig(filename)

    plt.close("all")
