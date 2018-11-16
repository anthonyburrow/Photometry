import numpy as np
from math import floor
from astropy import wcs
from astropy.io import fits
import csv

from .astrometry import GetAstrometryOffset


def ProcessMatch(cluster, date, app):
    """Matches targets between exposure times.

    Determines which targets have corresponding values between the short and
    long exposure data sets. The magnitudes used for each filter are determined
    by those with the lesser respective error. Targets with no matches are also
    included, as those typically refer to the brightest and darkest targets.

    Args:
        cluster (str): The cluster on which to perform matching.
        date (str): The date on which to perform matching.
        app (Application): The GUI application object that controls processing.

    Returns:
        list: 2-dimensional array consisting of X- and Y- image coordinates and
              the magnitudes and magnitude errors of each filter for every
              target.

    """
    print("Matching objects for " + date + "...\n")

    short_data = ProcessMatch_Filter(cluster, date, app, 'Short')
    long_data = ProcessMatch_Filter(cluster, date, app, 'Long')

    data = ProcessMatch_Exposure(cluster, date, app, short_data, long_data)

    data = ExtinctionCorrection(app, data)
    data = ApertureCorrection(date, cluster, data)

    GetRaDecs(cluster, date, data)

    # Output to file
    filename = 'output/' + cluster + '/' + date + \
               '/phot_' + app.phot_type + '.dat'
    with open(filename, 'w') as F:
        np.savetxt(F, [x[:-1] for x in data], fmt='%.3f')

    filename = 'output/' + cluster + '/' + date + \
               '/phot_exps_' + app.phot_type + '.dat'
    with open(filename, 'w') as F:
        for target in data:
            F.write('%8.3f' % target[0] + '    ')
            F.write('%8.3f' % target[1] + '    ')
            F.write('%6.3f' % target[2] + '    ')
            F.write('%6.3f' % target[4] + '    ')
            F.write('%6.3f' % target[6] + '    ')
            F.write('%6.3f' % target[8] + '    ')
            F.write(target[10])
            F.write("\n")

    return data


def alsRead(filename):
    """Reads PSF photometry files.

    PSF photometry is in the standard output .als format provided by the
    'allstar' task within IRAF's DAOPHOT package.  This file is located in
    `root/photometry/*date*/*cluster*/`.

    Args:
        filename (str): The input .als file (full path) to read.

    Returns:
        list: 2-dimensional array consisting of x- and y-image coordinates,
              magnitudes, and magnitude errors from the given file.

    """
    try:
        with open(filename) as F:
            file = F.readlines()[44:]
    except IOError:
        print("\nFile does not exist:\n" + filename)
        return

    data = []
    minimum = 0.
    maximum = 4096.
    for i in range(0, int(len(file) / 2.)):
        curr = file[2 * i].split()
        # Set image size limits (not needed if min/max equal image dimensions)
        if minimum < float(curr[1]) < maximum and \
           minimum < float(curr[2]) < maximum:
            combined = ' '.join([file[2 * i], file[2 * i + 1]])
            # Select values needed in data set: X, Y, mag, mag error
            combined = combined.split()
            selected = []
            selected.extend([float(combined[1]), float(combined[2]),
                             float(combined[3]), float(combined[4])])
            data.append(selected)

    return data


def magRead(filename):
    """Reads aperture photometry files.

    Aperture photometry is in the standard output .mag format provided by the
    phot' task within IRAF's DAOPHOT package.  This file is located in
    root/photometry/*date*/*cluster*/

    Args:
        filename (str): The input .mag file (full path) to read.

    Returns:
        list: 2-dimensional array consisting of x- and y-image coordinates,
              magnitudes, and magnitude errors.

    """
    try:
        with open(filename) as F:
            file = F.readlines()[75:]
    except IOError:
        print("\nFile does not exist:\n" + filename)
        return

    data = []

    for i in range(0, int(len(file) / 5.)):
        # Concatenate lines
        combined = ' '.join(file[5 * i: 5 * i + 5])
        # Select values needed in data set: X, Y, mag, mag error
        combined = combined.split()
        selected = []
        if combined[33] != 'INDEF' and combined[34] != 'INDEF':
            selected.extend((float(combined[7]), float(combined[8]),
                             float(combined[33]), float(combined[34])))
            data.append(selected)

    return data


def MatchTarget(tol, coords, data, indices=(0, 1)):
    """Primary matching system.

    This matching program matches a target with another target from a given data
    set that closest corresponds to it in x-y space.  Uses a pixel tolerance.
    The matched target is returned.

    Args:
        tol (float): The pixel tolerance for matching.
        coords (2-tuple): The x- and y-coordinates of the variable target
                          desired to be matched with.
        data (list): List of targets to be matched against.
        indices (2-tuple): Indices of `data` rows that contains x- and
                           y-coordinates desired to be matched with.  Defaults
                           to the first two indices.

    Returns:
        list: Returns a list corresponding to a row in `data` argument. Returns
              `None` if no match is found.

    """
    potentialMatches = []
    for target in data:
        if abs(coords[0] - target[indices[0]]) <= tol and \
           abs(coords[1] - target[indices[1]]) <= tol:
            potentialMatches.append(target)

    # Match with the closest target
    if potentialMatches:
        d = [np.sqrt((coords[0] - x[indices[0]])**2 +
                     (coords[1] - x[indices[1]])**2) for x in potentialMatches]
        d = np.array(d)
        match = potentialMatches[np.argmin(d)]
        return match

    return None


def ProcessMatch_Filter(cluster, date, app, exposure):
    """Matches targets between filters for each exposure time.

    The four filters used are B, V, R, and H-alpha.  Determines which targets
    have corresponding values within each filter data set.  Uses B filter data
    as a reference.

    Args:
        cluster (str): The cluster on which to perform matching.
        date (str): The date on which to perform matching.
        app (Application): The GUI application object that controls processing.
        exposure (str): Determines whether 'Short' or 'Long' exposure
                        times are used.

    Returns:
        list: 2-dimensional array consisting of x- and y-image coordinates and
              the magnitudes and magnitude errors for every target that has a
              corresponding value on each filter data set.

    """
    data = []

    if app.phot_type == 'psf':
        if exposure == 'Short':
            filenames = ['B1.als.1', 'V1.als.1', 'R1.als.1', 'H1.als.1']
        elif exposure == 'Long':
            filenames = ['B3.als.1', 'V3.als.1', 'R3.als.1', 'H3.als.1']

        path = 'photometry/' + cluster + '/' + date + '/'
        filenames = [path + x for x in filenames]

        B_data = alsRead(filenames[0])
        V_data = alsRead(filenames[1])
        R_data = alsRead(filenames[2])
        H_data = alsRead(filenames[3])

    elif app.phot_type == 'aperture':
        if exposure == 'Short':
            filenames = ['B1.mag.1', 'V1.mag.1', 'R1.mag.1', 'H1.mag.1']
        elif exposure == 'Long':
            filenames = ['B3.mag.1', 'V3.mag.1', 'R3.mag.1', 'H3.mag.1']

        path = 'photometry/' + cluster + '/' + date + '/'
        filenames = [path + x for x in filenames]

        B_data = magRead(filenames[0])
        V_data = magRead(filenames[1])
        R_data = magRead(filenames[2])
        H_data = magRead(filenames[3])

    # Get binning for observation night (assume same bin per night)
    filename = 'photometry/' + cluster + '/' + date + '/B1.fits'
    F = fits.getheader(filename)
    binning = F['XBINNING']

    # Specify any coordinate offsets left to be made (from the B image, which
    # is the reference)
    if exposure == 'Short':
        V_coo_offset = GetAstrometryOffset(cluster, date, baseDate=date,
                                           image='V1', baseImage='B1') / binning
        R_coo_offset = GetAstrometryOffset(cluster, date, baseDate=date,
                                           image='R1', baseImage='B1') / binning
        H_coo_offset = GetAstrometryOffset(cluster, date, baseDate=date,
                                           image='H1', baseImage='B1') / binning
    if exposure == 'Long':
        V_coo_offset = GetAstrometryOffset(cluster, date, baseDate=date,
                                           image='V3', baseImage='B3') / binning
        R_coo_offset = GetAstrometryOffset(cluster, date, baseDate=date,
                                           image='R3', baseImage='B3') / binning
        H_coo_offset = GetAstrometryOffset(cluster, date, baseDate=date,
                                           image='H3', baseImage='B3') / binning

    # Apply coordinate corrections
    for target in V_data:
        target[0] += V_coo_offset[0]
        target[1] += V_coo_offset[1]
    for target in R_data:
        target[0] += R_coo_offset[0]
        target[1] += R_coo_offset[1]
    for target in H_data:
        target[0] += H_coo_offset[0]
        target[1] += H_coo_offset[1]

    # Match stars between filters
    data = []
    for target in B_data:
        v_match = MatchTarget(app.cooTol / binning, [target[0], target[1]],
                              V_data)
        if v_match is None:
            continue
        r_match = MatchTarget(app.cooTol / binning, [target[0], target[1]],
                              R_data)
        if r_match is None:
            continue
        h_match = MatchTarget(app.cooTol / binning, [target[0], target[1]],
                              H_data)
        if h_match is None:
            continue
        target.extend(v_match + r_match + h_match)
        data.append(target)

        V_data.remove(v_match)   #
        R_data.remove(r_match)   # Stop from matching again
        H_data.remove(h_match)   #

    if exposure == 'Short':
        for target in data:
            target.append('s')
    elif exposure == 'Long':
        for target in data:
            target.append('l')

    print("    " + exposure + " matched: " + str(len(data)))

    return data


def ProcessMatch_Exposure(cluster, date, app, short_data, long_data):
    """Matches targets between exposure times.

    Uses matching program to find corresponding targets between long- and
    short-exposure data sets and uses the photometry for the target with
    less error or that of the unsaturated target.

    Args:
        cluster (str): The cluster on which to perform matching.
        date (str): The date on which to perform matching.
        app (Application): The GUI application object that controls processing.
        short_data (list): The data set of short-exposure targets.
        long_data (list): The data set of long-exposure targets.

    Returns:
        list: 2-dimensional array consisting of x- and y-image coordinates and
              the magnitudes and magnitude errors for every target after
              completed matching.

    """
    # Get binning for observation night (assume same bin per night)
    filename = 'photometry/' + cluster + '/' + date + '/B1.fits'
    F = fits.getheader(filename)
    binning = F['XBINNING']

    # Get coordinate offset between B1 and B3 (short and long exp)
    coord_offset = GetAstrometryOffset(cluster, date, baseDate=date, image='B3',
                                       baseImage='B1') / binning

    # Match between short and long exposures and use values from that with the
    # lowest error
    print("\n  Matching objects between long and short exposures...")

    data = short_data + long_data
    matches = []
    count = 0

    for target in short_data:
        coords = [target[0] - coord_offset[0], target[1] - coord_offset[1]]
        match = MatchTarget(app.cooTol / binning, coords, long_data)

        if match is None:
            continue

        exps = ""
        phot_used = [target[0], target[1]]   # Use B short coordinates

        for i, fil in enumerate(['B', 'V', 'R', 'H']):
            x_i = 4 * i
            y_i = x_i + 1
            pho_i = x_i + 2
            err_i = x_i + 3
            x = floor(match[x_i])
            y = floor(match[y_i])
            if (target[err_i] >= match[err_i]) and \
               SaturationCheck(cluster, date, fil, x, y):
                phot_used.extend([match[pho_i], match[err_i]])
                exps += "l"
            else:
                phot_used.extend([target[pho_i], target[err_i]])
                exps += "s"

        phot_used.append(exps)

        data.remove(target)
        data.remove(match)
        long_data.remove(match)  # Stop from matching again
        matches.append(phot_used)

        count += 1

    data = [[x[0], x[1], x[2], x[3], x[6], x[7],
             x[10], x[11], x[14], x[15], x[16]] for x in data]
    data.extend(matches)

    print("\n    Matched between exposures: " + str(count))
    print("    Short only: " + str(len(short_data) - count))
    print("    Long only: " + str(len(long_data) - count))
    print("    Total: " + str(len(data)))

    return data


def SaturationCheck(cluster, date, fil, x, y):
    """Determines if a point on an image is saturated.

    Args:
        cluster (str): The cluster the to which the image corresponds.
        date (str): The date the to which the image corresponds.
        fil (str): String corresponding to the filter (B, V, R, or H) of
                   the image.
        x (float): x-coordinate on the image.
        y (float): y-coordinate on the image.

    Returns:
        bool: Returns True if point in image is not saturated, and False
              if it is saturated.

    """
    saturation = 60000

    filename = 'photometry/' + cluster + '/' + date + '/' + fil + '3.fits'
    with fits.open(filename) as hdu:
        data = hdu[0].data
        val = data[y, x]
        # Return true if not saturated
        return val < saturation


def GetRaDecs(cluster, date, data):
    """Creates a list of RA/Decs corresponding to list of matched targets.

    Args:
        cluster (str): The cluster from which the corresponding `wcs` is
                       located.
        date (str): The date from which the corresponding `wcs` is located.
        data (list): Data set from which the RA/Decs should be extracted

    """
    w = wcs.WCS('photometry/' + cluster + '/' + date + '/B1_wcs.fits')
    coords = []
    for target in data:
        radec = w.all_pix2world(target[0], target[1], 0)
        ra = float(radec[0])
        dec = float(radec[1])
        coords.append([ra, dec])

    filename = 'photometry/' + cluster + '/' + date + '/phot_radec.csv'
    with open(filename, 'w') as F:
        writer = csv.writer(F)
        writer.writerows(coords)


def ExtinctionCorrection(app, data):
    """Applies extinction correction.

    Args:
        app (Application): The GUI application object that controls processing.
        data (list): Data set on which to apply corrections.

    Returns:
        list: Data set with corrected photometry.

    """
    for target in data:
        target[2] -= app.A_b
        target[4] -= app.A_v
        target[6] -= app.A_r

    return data


def ApertureCorrection(date, cluster, data):
    """Applies extinction correction.

    Args:
        date (str): Date corresponding to specific aperture corrections.
        cluster (str): Cluster corresponding to specific aperture corrections.
        data (list): Data set on which to apply corrections.

    Returns:
        list: Data set with corrected photometry.

    """
    try:
        filename = 'standards/' + date + '/' + cluster + \
                   '_aperture_corrections.dat'
        corrections = np.loadtxt(filename)
    except IOError:
        print("  \nAperture corrections not applied.")
        return data

    for target in data:
        target[2] += corrections[0]
        target[4] += corrections[1]
        target[6] += corrections[2]
    print("  \nAperture corrections applied.\n")

    return data
