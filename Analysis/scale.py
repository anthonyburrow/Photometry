import numpy as np
from astropy.io import fits

import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator, FormatStrFormatter

from os import remove
import os.path

from .observations import ListDates
from .astrometry import GetAstrometryOffset


def SetData(cluster, date, app):
    """Determines which targets are not within a given tolerance of another.

    Given an input data set, this outputs those targets which are not within a
    specified spacial tolerance in pixels.

    Args:
            date:

    Returns:
            2-dimensional array consisting of X- and Y- image coordinates, magnitudes,
            and magnitude errors.

    """
    binning = Binning(cluster, date)

    # Read low error data, or else read normal data
    try:
        filename = 'output/' + cluster + '/' + date + '/phot_' + app.phot_type + '_lowError.dat'
        data = np.loadtxt(filename).tolist()
    except IOError:
        try:
            filename = 'output/' + cluster + '/' + date + '/phot_' + app.phot_type + '.dat'
            data = np.loadtxt(filename).tolist()
        except IOError:
            print("\nFile does not exist:\n" + filename)
            return

    # Get rid of stars with other stars next to them
    for target in reversed(data):
        for otherTarget in data:
            r = binning * np.sqrt((target[0] - otherTarget[0])**2 + (target[1] - otherTarget[1])**2)
            if r <= app.cooTol and target != otherTarget:
                data.remove(target)
                break

    # Get rid of stars that are outliers
    try:
        filename = 'output/' + cluster + '/' + date + '/beList_' + app.phot_type + '.dat'
        filtered_data = np.loadtxt(filename).tolist()
        for target in reversed(data):
            if target in filtered_data:
                data.remove(target)
    except IOError:
        print("  Note: Outliers were not removed from scale sample because 'beList_" + app.phot_type + ".dat' does not exist.")

    return data


def Binning(cluster, date):
    """Determines the bin size used during the observation.

    Reads the B1.fits image for the specified date of observation and extracts
    the binning data.  This uses x-axis binning, however x- and y- axis binning
    are typically equal.

    Args:
            date: The date of observation.

    Returns:
            The bin size as an integer.

    """
    filename = 'photometry/' + cluster + '/' + date + '/B1.fits'
    try:
        F = fits.getheader(filename)
        binning = F['XBINNING']
    except IOError:
        print("\nFile does not exist:\n" + filename)
        binning = 1

    return binning


def ProcessScale(cluster, date, app, baseDate):
    """Scales the given set of data.

    Corresponds each target in the data set to others in the reference data set and
    averages magnitude differences between each set of targets to create an averaged
    magnitude offset for each filter.  Does not use either targets within a certain
    x-y distance of any others or targets that are considered outliers.

    """
    # Reference date process
    if date == baseDate:
        # Set documented scale offsets to zero
        offsets = [[0, 0], [0, 0], [0, 0], [0, 0]]
        with open('output/' + cluster + '/' + date + '/magScales.dat', 'w') as F:
            np.savetxt(F, offsets, fmt="%.3f")

        # Remove any unneeded/extra files (usually after re-scaling)
        path = 'output/' + cluster + '/' + date + '/'
        files = [
            path + "plots/num_vs_Bmag_diffs_" + app.phot_type + ".png",
            path + "plots/num_vs_Vmag_diffs_" + app.phot_type + ".png",
            path + "plots/num_vs_Rmag_diffs_" + app.phot_type + ".png",
            path + "plots/num_vs_Hmag_diffs_" + app.phot_type + ".png"
        ]
        for file in files:
            if os.path.isfile(file):
                remove(file)

        return

    print("\nScaling data for " + cluster + " on " + date + "...\n")

    # Read data
    print("  Creating reference sample...")
    baseBinning = Binning(cluster, baseDate)
    baseData = SetData(cluster, baseDate, app)

    print("  Creating variable sample...")
    binning = Binning(cluster, date)
    data = SetData(cluster, date, app)

    print("\n")

    # Get x- and y-offsets
    xOffset, yOffset = GetAstrometryOffset(cluster, date, baseDate)

    # Find corresponding targets
    B_diff = []
    V_diff = []
    R_diff = []
    H_diff = []

    for target in data:
        for baseTarget in baseData:
            if abs(baseBinning * baseTarget[0] - (binning * target[0] + xOffset)) <= app.cooTol and \
               abs(baseBinning * baseTarget[1] - (binning * target[1] + yOffset)) <= app.cooTol:
                B_diff.append(baseTarget[2] - target[2])
                V_diff.append(baseTarget[4] - target[4])
                R_diff.append(baseTarget[6] - target[6])
                H_diff.append(baseTarget[8] - target[8])

    offsets = GetOffsets([B_diff, V_diff, R_diff, H_diff])

    # Print scale information
    print("\n  Scaled with ", len(B_diff), " stars:")
    print("  B offset = " + "%.3f" % offsets[0][0] + " +/- " + "%.3f" % offsets[0][1])
    print("  V offset = " + "%.3f" % offsets[1][0] + " +/- " + "%.3f" % offsets[1][1])
    print("  R offset = " + "%.3f" % offsets[2][0] + " +/- " + "%.3f" % offsets[2][1])
    print("  H offset = " + "%.3f" % offsets[3][0] + " +/- " + "%.3f" % offsets[3][1])

    # Output to file
    with open('output/' + cluster + '/' + date + '/magScales.dat', 'w') as F:
        np.savetxt(F, offsets, fmt='%.3f')


def GetOffsets(diffs):
    offsets = []
    for diff in diffs:
        if diff:
            offset = []
            offset.append(np.mean(diff))
            offset.append(np.std(diff))

            # Recalculate using only the mag differences within 3-sigma
            for target in reversed(diff):
                if not offset[0] - 3 * offset[1] < target < offset[0] + 3 * offset[1]:
                    diff.remove(target)
            offset[0] = np.mean(diff)
            offset[1] = np.std(diff)
        else:
            offset = [0, 0]

        offsets.append(offset)

        titles = ['B', 'V', 'R', 'H-alpha']
        # self.num_vs_mag_hist(cluster, date, diff, offset[0], offset[1], titles[diffs.index(diff)])

    return offsets


def num_vs_mag_hist(cluster, date, app, x, mean, std, filter):
    plt.figure(figsize=(12, 9))

    # Plot main data
    plt.hist(x, bins=20, range=(mean - 3 * std, mean + 3 * std), color='#3f3f3f')

    # plt.title(filter + " Magnitude Scaling Differences")
    plt.xlabel('Magnitude Difference', fontsize=36)
    plt.ylabel('Frequency', fontsize=36)

    plt.axes().xaxis.set_major_locator(MultipleLocator(0.1))
    plt.axes().xaxis.set_major_formatter(FormatStrFormatter('%.2f'))
    plt.axes().xaxis.set_minor_locator(MultipleLocator(0.025))

    # plt.axes().yaxis.set_major_locator(MultipleLocator(10))
    # plt.axes().yaxis.set_major_formatter(FormatStrFormatter('%d'))
    # plt.axes().yaxis.set_minor_locator(MultipleLocator(2.5))

    plt.axes().tick_params('both', length=12, width=4, which='major', top=True, right=True, direction='in', pad=6, labelsize=30)
    plt.axes().tick_params('both', length=8, width=3, which='minor', top=True, right=True, direction='in')

    plt.axes().spines['top'].set_linewidth(4)
    plt.axes().spines['right'].set_linewidth(4)
    plt.axes().spines['bottom'].set_linewidth(4)
    plt.axes().spines['left'].set_linewidth(4)

    ymin, ymax = plt.ylim()
    plt.vlines(mean + std, 0, ymax, colors='#ff5151', linestyles='dashed', label='Standard Error')
    plt.vlines(mean - std, 0, ymax, colors='#ff5151', linestyles='dashed')

    # plt.legend(fontsize=28)
    plt.tight_layout()

    # Output
    filename = 'output/' + cluster + '/' + date + '/plots/num_vs_' + filter[0] + 'mag_diffs_' + app.phot_type + '.png'
    plt.savefig(filename)

    plt.clf()


def Rescale(cluster, app):
    print("\nFinding optimal observation date with which to re-scale data...")

    # Look for "brightest" night
    check = []
    dates = ListDates(cluster)
    for date in dates:
        filename = 'output/' + cluster + '/' + date + '/magScales.dat'
        scales = np.loadtxt(filename)
        check.append(scales[1][0])   # check against V mag just because

    brightestIndex = check.index(max(check))
    newBaseDate = ListDates(cluster)[brightestIndex]

    # Rerun scale process with new reference date
    print("  Re-scaling data with reference date " + newBaseDate)
    for date in dates:
        ProcessScale(cluster, date, app, newBaseDate)
        ApplyScale(cluster, date, app)


def ApplyScale(cluster, date, app):
    # Implement scale offsets
    filename = 'output/' + cluster + '/' + date + '/phot_' + app.phot_type + '.dat'
    orig_data = np.loadtxt(filename)

    filename = 'output/' + cluster + '/' + date + '/magScales.dat'
    scales = np.loadtxt(filename)

    for target in orig_data:
        target[2] += scales[0][0]
        target[4] += scales[1][0]
        target[6] += scales[2][0]
        target[8] += scales[3][0]

        # Decided against adding scaling error contribution
        # target[3] = np.sqrt(target[3]**2 + scales[0][1]**2)
        # target[5] = np.sqrt(target[5]**2 + scales[1][1]**2)
        # target[7] = np.sqrt(target[7]**2 + scales[2][1]**2)
        # target[9] = np.sqrt(target[9]**2 + scales[3][1]**2)

    # Write to file
    filename = 'output/' + cluster + '/' + date + '/phot_scaled_' + app.phot_type + '.dat'
    with open(filename, 'w') as F:
        np.savetxt(F, orig_data, fmt="%.3f")
