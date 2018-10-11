import numpy as np
from astropy import wcs
from astropy.io import fits
import csv

from .astrometry import GetAstrometryOffset


def ProcessMatch(cluster, date, app):
    """Matches targets between exposure times.

    Determines which targets have corresponding values between the short and
    long exposure data sets.  The magnitudes used for each filter are determined
    by those with the lesser respective error.  Targets with no matches are also
    included, as those typically refer to the brightest and darkest targets.

    Returns:
            2-dimensional array consisting of X- and Y- image coordinates and the
            magnitudes and magnitude errors of each filter for every target.

    """
    print("Matching objects for " + date + "...\n")

    short_data = ProcessMatch_Filter(cluster, date, app, 'Short')
    long_data = ProcessMatch_Filter(cluster, date, app, 'Long')

    data = ProcessMatch_Exposure(cluster, date, app, short_data, long_data)

    data = ExtinctionCorrection(app, data)
    data = ApertureCorrection(date, cluster, data)

    GetRaDecs(cluster, date, data)

    # Output to file
    filename = 'output/' + cluster + '/' + date + '/phot_' + app.phot_type + '.dat'
    with open(filename, 'w') as F:
        np.savetxt(F, [x[:-1] for x in data], fmt='%.3f')

    filename = 'output/' + cluster + '/' + date + '/phot_exps_' + app.phot_type + '.dat'
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

    PSF photometry is in the standard output .als format provided by the 'allstar'
    task within IRAF's DAOPHOT package.  This file is located in
    root/photometry/*date*/*cluster*/

    Args:
            file (string): The input .als file to read.

    Returns:
            2-dimensional array consisting of X- and Y- image coordinates, magnitudes,
            and magnitude errors.

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
            selected.extend([float(combined[1]), float(combined[2]), float(combined[3]), float(combined[4])])
            data.append(selected)

    return data


def magRead(filename):
    """Reads aperture photometry files.

    Aperture photometry is in the standard output .mag format provided by the 'phot'
    task within IRAF's DAOPHOT package.  This file is located in
    root/photometry/*date*/*cluster*/

    Args:
            file (string): The input .mag file to read.

    Returns:
            2-dimensional array consisting of X- and Y- image coordinates, magnitudes,
            and magnitude errors.

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
            selected.extend((float(combined[7]), float(combined[8]), float(combined[33]), float(combined[34])))
            data.append(selected)

    return data


def MatchTarget(tol, coords, data, indices=[0, 1]):
    potentialMatches = []
    for target in data:
        if abs(coords[0] - target[indices[0]]) <= tol and abs(coords[1] - target[indices[1]]) <= tol:
            potentialMatches.append(target)

    # Match with the closest target
    if potentialMatches:
        d = [np.sqrt(coords[0] - x[indices[0]])**2 + (coords[1] - x[indices[1]])**2 for x in potentialMatches]
        d = np.array(d)
        match = potentialMatches[np.argmin(d)]
        return match

    return None


def ProcessMatch_Filter(cluster, date, app, exposure):
    """Matches targets between filters for each exposure time.

    The four filters used are B, V, R, and H-alpha.  Determines which targets have
    corresponding values within each filter data set.  Uses B filter data as a
    reference.

    Args:
            exposure (string): Determines whether "Short" or "Long" exposure times are used.

    Returns:
            2-dimensional array consisting of X- and Y- image coordinates and the
            magnitudes and magnitude errors for every target that has a corresponding
            value on each filter data set.

    """
    # Create data sets for each filter
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

    # Specify any coordinate offsets left to be made (from the B image, which is the reference)
    if exposure == 'Short':
        V_coo_offset = GetAstrometryOffset(cluster, date, baseDate=date, image='V1', baseImage='B1') / binning
        R_coo_offset = GetAstrometryOffset(cluster, date, baseDate=date, image='R1', baseImage='B1') / binning
        H_coo_offset = GetAstrometryOffset(cluster, date, baseDate=date, image='H1', baseImage='B1') / binning
    if exposure == 'Long':
        V_coo_offset = GetAstrometryOffset(cluster, date, baseDate=date, image='V3', baseImage='B3') / binning
        R_coo_offset = GetAstrometryOffset(cluster, date, baseDate=date, image='R3', baseImage='B3') / binning
        H_coo_offset = GetAstrometryOffset(cluster, date, baseDate=date, image='H3', baseImage='B3') / binning

    print(exposure + ' V offset: [%s, %s]' % (V_coo_offset[0], V_coo_offset[1]))
    print(exposure + ' R offset: [%s, %s]' % (R_coo_offset[0], R_coo_offset[1]))
    print(exposure + ' H offset: [%s, %s]' % (H_coo_offset[0], H_coo_offset[1]))

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
        v_match = MatchTarget(app.cooTol / binning, [target[0], target[1]], V_data)
        if v_match is None:
            continue
        r_match = MatchTarget(app.cooTol / binning, [target[0], target[1]], R_data)
        if r_match is None:
            continue
        h_match = MatchTarget(app.cooTol / binning, [target[0], target[1]], H_data)
        if h_match is None:
            continue
        target.extend([v_match[2], v_match[3],
                       r_match[2], r_match[3],
                       h_match[2], h_match[3]])
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
    # Get binning for observation night (assume same bin per night)
    filename = 'photometry/' + cluster + '/' + date + '/B1.fits'
    F = fits.getheader(filename)
    binning = F['XBINNING']

    # Apply coordinate offset between B1 and B3
    coord_offset = GetAstrometryOffset(cluster, date, baseDate=date, image='B3', baseImage='B1') / binning
    for target in long_data:
        target[0] += coord_offset[0]
        target[1] += coord_offset[1]

    print('\nLong offset: [%s, %s]' % (coord_offset[0], coord_offset[1]))

    # Match between short and long exposures and use values from that with the lowest error
    print("\n  Matching objects between long and short exposures...")

    data = short_data + long_data
    count = 0

    for target in short_data:
        match = MatchTarget(app.cooTol / binning, [target[0], target[1]], long_data)

        if match is None:
            continue

        exps = ""
        phot_used = [target[0], target[1]]

        for n in [0, 2, 4, 6]:
            if (target[3 + n] >= match[3 + n]) and (abs(target[2 + n] - match[2 + n]) < app.magTol):
                phot_used.extend((match[2 + n], match[3 + n]))
                exps += "l"
            else:
                phot_used.extend((target[2 + n], target[3 + n]))
                exps += "s"

        phot_used.append(exps)

        data.remove(target)
        data.remove(match)
        long_data.remove(match)  # Stop from matching again
        data.append(phot_used)

        count += 1

    print("\n    Matched between exposures: " + str(count))
    print("    Short only: " + str(len(short_data) - count))
    print("    Long only: " + str(len(long_data) - count))
    print("    Total: " + str(len(data)))

    return data


def GetRaDecs(cluster, date, data):
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
    for target in data:
        target[2] -= app.A_b
        target[4] -= app.A_v
        target[6] -= app.A_r

    return data


def ApertureCorrection(date, cluster, data):
    try:
        filename = 'standards/' + date + '/' + cluster + '_aperture_corrections.dat'
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
