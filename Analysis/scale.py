import numpy as np

from .io import Binning
from .observations import ListDates
from .astrometry import GetAstrometryOffset
from .low_error import ProcessLowError


def SetData(cluster, date, app):
    """Determines which targets are not within a given tolerance of another.

    Given an input data set, this outputs those targets which are not within a
    specified spacial tolerance in pixels, and which are not Be candidates
    (photometric outliers).

    Args:
        cluster (str): Cluster from which data is retrieved.
        date (str): Date from which data is retrieved.
        app (Application): The GUI application object that controls processing.

    Returns:
        list: List of data suitable for sampling.

    """
    binning = Binning(cluster, date)

    # Use low-error data
    filename = 'output/%s/%s/phot.dat' % (cluster, date)
    data = np.loadtxt(filename)
    data = ProcessLowError(cluster, date, data).tolist()

    # Get rid of stars with other stars next to them
    for target in reversed(data):
        for otherTarget in data:
            r = binning * np.sqrt((target[0] - otherTarget[0])**2 +
                                  (target[1] - otherTarget[1])**2)
            if r <= app.cooTol and target != otherTarget:
                data.remove(target)
                break

    # Get rid of stars that are outliers
    try:
        filename = 'output/%s/%s/belist.dat' % (cluster, date)
        filtered_data = np.loadtxt(filename).tolist()
        for target in reversed(data):
            if target in filtered_data:
                data.remove(target)
    except IOError:
        print("  Note: Outliers were not removed from scale sample because \
               'beList.dat' does not exist.")

    return data


def ProcessScale(cluster, date, app, baseDate):
    """Controls the full process of scaling data.

    Corresponds each target in the data set to others in the reference data set
    and averages magnitude differences between each set of targets to create an
    averaged magnitude offset for each filter.  Does not use either targets
    within a certain x-y distance of any others or targets that are considered
    outliers.

    Args:
        cluster (str): Cluster from which data is retrieved.
        date (str): Date from which data is retrieved.
        app (Application): The GUI application object that controls processing.
        baseDate (str): Reference date against which data is scaled.

    """
    # Reference date process
    if date == baseDate:
        # Set documented scale offsets to zero
        offsets = [[0, 0], [0, 0], [0, 0], [0, 0]]
        filename = 'output/%s/%s/mag_scales.dat' % (cluster, date)
        with open(filename, 'w') as F:
            np.savetxt(F, offsets, fmt="%.3f")

        return

    # Read data
    print("  Creating reference sample...")
    baseBinning = Binning(cluster, baseDate)
    baseData = SetData(cluster, baseDate, app)

    print("  Creating variable sample...\n")
    binning = Binning(cluster, date)
    data = SetData(cluster, date, app)

    # Get x- and y-offsets
    xOffset, yOffset = GetAstrometryOffset(cluster, date, baseDate)

    # Find corresponding targets
    B_diff = []
    V_diff = []
    R_diff = []
    H_diff = []

    for target in data:
        for baseTarget in baseData:
            if abs(baseBinning * baseTarget[0] -
                   (binning * target[0] + xOffset)) <= app.cooTol and \
               abs(baseBinning * baseTarget[1] -
                   (binning * target[1] + yOffset)) <= app.cooTol:
                B_diff.append(baseTarget[2] - target[2])
                V_diff.append(baseTarget[4] - target[4])
                R_diff.append(baseTarget[6] - target[6])
                H_diff.append(baseTarget[8] - target[8])

    offsets = GetOffsets([B_diff, V_diff, R_diff, H_diff])

    # Output
    t = '  Scaled with %d stars:\n' % len(B_diff) + \
        '    B offset = %.3f +/- %.3f\n' % (offsets[0][0], offsets[0][1]) + \
        '    V offset = %.3f +/- %.3f\n' % (offsets[1][0], offsets[1][1]) + \
        '    R offset = %.3f +/- %.3f\n' % (offsets[2][0], offsets[2][1]) + \
        '    H offset = %.3f +/- %.3f\n' % (offsets[3][0], offsets[3][1])

    print(t)

    filename = 'output/%s/%s/mag_scales.dat' % (cluster, date)
    with open(filename, 'w') as F:
        np.savetxt(F, offsets, fmt='%.3f')


def GetOffsets(diffs):
    """Calculates mean photometric offsets given a set of individual photometry
    offsets.

    Calculates the mean of a given set, then recalculates based on the
    calculated 3-sigma.

    Args:
        diffs (list): List of sets of individual photometric offsets for each
                      filter.

    Returns:
        list: Set of mean photometric offsets for each filter.

    """
    offsets = []
    for diff in diffs:
        if diff:
            offset = []
            offset.append(np.mean(diff))
            offset.append(np.std(diff))

            # Recalculate using only the mag differences within 3-sigma
            for target in reversed(diff):
                if not abs(target - offset[0]) < 3 * offset[1]:
                    diff.remove(target)
            offset[0] = np.mean(diff)
            offset[1] = np.std(diff)
        else:
            offset = [0, 0]

        offsets.append(offset)

    return offsets


def Rescale(cluster, app):
    """Rescales the data based on a different reference date.

    Reference date is found by using the date with the highest photometry
    offsets, corresponding to the brightest image, meaning less atmospheric
    intrusion.

    Args:
        cluster (str): Cluster from which data is retrieved.
        app (Application): The GUI application object that controls processing.

    """
    # Look for "brightest" night
    check = []
    dates = ListDates(cluster)
    for date in dates:
        filename = 'output/%s/%s/mag_scales.dat' % (cluster, date)
        scales = np.loadtxt(filename)
        check.append(scales[1][0])   # check against V mag

    brightestIndex = check.index(max(check))
    newBaseDate = ListDates(cluster)[brightestIndex]

    # Rerun scale process with new reference date
    for date in dates:
        print("Scaling data with final scaling for %s on %s using " %
              (cluster, date) + "reference %s...\n" % newBaseDate)
        ProcessScale(cluster, date, app, newBaseDate)
        ApplyScale(cluster, date)


def ApplyScale(cluster, date):
    """Applies photometric offset corrections to original data.

    Args:
        cluster (str): Cluster from which data is retrieved.
        date (str): Date from which data is retrieved.

    """
    filename = 'output/%s/%s/mag_scales.dat' % (cluster, date)
    scales = np.loadtxt(filename)

    for file in ('', '_accepted', '_rejected'):
        filename = 'output/%s/%s/phot%s.dat' % (cluster, date, file)
        orig_data = np.loadtxt(filename)

        for target in orig_data:
            target[2] += scales[0][0]
            target[4] += scales[1][0]
            target[6] += scales[2][0]
            target[8] += scales[3][0]

        # Write to file
        filename = 'output/%s/%s/phot_scaled%s.dat' % (cluster, date, file)
        with open(filename, 'w') as F:
            np.savetxt(F, orig_data, fmt="%.3f")
