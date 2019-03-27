import numpy as np
from math import floor
from astropy.io import fits
import csv
import os.path

from .io import Binning, GetWCS, magRead
from .astrometry import GetAstrometryOffset

import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec


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
    print("  Matching objects between long and short exposures...\n")
    short_data = ProcessMatch_Filter(cluster, date, app, 'Short')
    long_data = ProcessMatch_Filter(cluster, date, app, 'Long')

    print("\n  Matching objects between long and short exposures...\n")
    data = ProcessMatch_Exposure(cluster, date, app, short_data, long_data)

    print("\n  Applying extinction correction...")
    data = ExtinctionCorrection(cluster, data)

    print("  Applying aperture correction...\n")
    data = ApertureCorrection(date, cluster, data)

    GetRaDecs(cluster, date, data)

    # Output to file
    filename = 'output/%s/%s/phot.dat' % (cluster, date)
    with open(filename, 'w') as F:
        np.savetxt(F, [x[:-1] for x in data], fmt='%.3f')

    filename = 'output/%s/%s/phot_exps.dat' % (cluster, date)
    with open(filename, 'w') as F:
        for target in data:
            t = '%8.3f' % target[0] + '\t' + \
                '%8.3f' % target[1] + '\t' + \
                '%6.3f' % target[2] + '\t' + \
                '%6.3f' % target[5] + '\t' + \
                '%6.3f' % target[8] + '\t' + \
                '%6.3f' % target[11] + '\t' + \
                target[14] + \
                '\n'
            F.write(t)

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

    if exposure == 'Short':
        filenames = ('B1.mag.1', 'V1.mag.1', 'R1.mag.1', 'H1.mag.1')
    elif exposure == 'Long':
        filenames = ('B3.mag.1', 'V3.mag.1', 'R3.mag.1', 'H3.mag.1')

    path = 'photometry/%s/%s/' % (cluster, date)
    filenames = [path + x for x in filenames]

    B_data = magRead(filenames[0])
    V_data = magRead(filenames[1])
    R_data = magRead(filenames[2])
    H_data = magRead(filenames[3])

    # Get binning for observation night (assume same bin per night)
    binning = Binning(cluster, date)

    # Specify any coordinate offsets left to be made (from the B image, which
    # is the reference)
    if exposure == 'Short':
        V_coo_offset = GetAstrometryOffset(cluster, date, baseDate=date,
                                           image='V1', baseImage='B1')
        R_coo_offset = GetAstrometryOffset(cluster, date, baseDate=date,
                                           image='R1', baseImage='B1')
        H_coo_offset = GetAstrometryOffset(cluster, date, baseDate=date,
                                           image='H1', baseImage='B1')
    if exposure == 'Long':
        V_coo_offset = GetAstrometryOffset(cluster, date, baseDate=date,
                                           image='V3', baseImage='B3')
        R_coo_offset = GetAstrometryOffset(cluster, date, baseDate=date,
                                           image='R3', baseImage='B3')
        H_coo_offset = GetAstrometryOffset(cluster, date, baseDate=date,
                                           image='H3', baseImage='B3')

    # Match stars between filters
    data = []
    for target in B_data:
        x = target[0] - (V_coo_offset[0] / binning)
        y = target[1] - (V_coo_offset[1] / binning)
        v_match = MatchTarget(app.cooTol / binning, (x, y), V_data)
        if v_match is None:
            continue

        x = target[0] - (R_coo_offset[0] / binning)
        y = target[1] - (R_coo_offset[1] / binning)
        r_match = MatchTarget(app.cooTol / binning, (x, y), R_data)
        if r_match is None:
            continue

        x = target[0] - (H_coo_offset[0] / binning)
        y = target[1] - (H_coo_offset[1] / binning)
        h_match = MatchTarget(app.cooTol / binning, (x, y), H_data)
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

    print("    %s matched: %d" % (exposure, len(data)))

    return data


def get_counts(mag, exptime):
    ZP = 21
    flux = exptime * 10**((mag - ZP) / -2.5)
    return flux


def get_exptime(cluster, date, fil, exposure):
    if exposure == 'Short':
        ex = '1'
    elif exposure == 'Long':
        ex = '3'

    filename = 'photometry/%s/%s/%s%s.fits' % (cluster, date, fil, ex)
    try:
        F = fits.getheader(filename)
        exptime = F['EXPTIME']
    except IOError:
        print("\nFile does not exist:\n%s" % filename)
        return

    return exptime


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
    binning = Binning(cluster, date)

    # Get coordinate offset between B1 and B3 (short and long exp)
    coord_offset = GetAstrometryOffset(cluster, date, baseDate=date, image='B3',
                                       baseImage='B1')
    coord_offset = (coord_offset[0] / binning, coord_offset[1] / binning)

    # Match between short and long exposures and use values from that with the
    # lowest error
    data = short_data + long_data
    matches = []
    count = 0

    plt.style.use('researchpaper')
    fig = plt.figure(dpi=200)
    fig.set_size_inches(20, 9)

    gs = GridSpec(1, 2)

    ax1 = fig.add_subplot(gs[0, 0])
    ax2 = fig.add_subplot(gs[0, 1])

    plot_colors = {
        'B': ('#028ecc', '#5ec3f2'),   # blue
        'V': ('#b933ea', '#c877e5'),   # purple
        'R': ('#ea3535', '#ff8787'),   # red
        'H': ('#ea8f34', '#ffbf87'),   # orange
    }

    for target in short_data:
        coords = (target[0] - coord_offset[0], target[1] - coord_offset[1])
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

            s_mag = target[pho_i]
            l_mag = match[pho_i]
            s_err = target[err_i]
            l_err = match[err_i]

            s_exptime = get_exptime(cluster, date, fil, 'Short')
            l_exptime = get_exptime(cluster, date, fil, 'Long')

            s_counts = get_counts(s_mag, s_exptime)
            l_counts = get_counts(l_mag, l_exptime)

            l_rate = l_counts / l_exptime
            s_rate = s_counts / s_exptime
            x1 = l_rate   # Short
            y1 = s_rate / l_rate   # Long - Short

            x2 = l_mag
            y2 = l_mag - s_mag
            x2err = s_err
            y2err = np.sqrt(s_err**2 + l_err**2)

            saturated = SaturationCheck(cluster, date, fil, match[x_i], match[y_i])

            mk = 'o'
            if saturated:
                mk = 'x'

            ax1.plot(x1, y1, mk, color=plot_colors[fil][0],
                     markersize=11, alpha=0.7)

            ax2.plot(x2, y2, mk, color=plot_colors[fil][0],
                     markersize=11, alpha=0.7)
            ax2.errorbar(x2, y2, xerr=x2err, yerr=y2err, fmt='none',
                         ecolor=plot_colors[fil][1], elinewidth=7, alpha=0.7)

            # Short -> Long correction
            s_err_high = s_err
            s_err_low = s_err
            if not saturated:
                ls_offset = l_mag - s_mag
                if ls_offset > 0:
                    s_err_high += ls_offset
                elif ls_offset < 0:
                    s_err_low += ls_offset

            if (s_err >= l_err) and not saturated:
                phot_used.extend([l_mag, l_err, l_err])
                exps += "l"
            else:
                phot_used.extend([s_mag, s_err_low, s_err_high])
                exps += "s"

        phot_used.append(exps)

        data.remove(target)
        data.remove(match)
        long_data.remove(match)  # Stop from matching again
        matches.append(phot_used)

        count += 1

    for target in data:
        if target in long_data:
            target[0] += coord_offset[0]
            target[1] += coord_offset[1]

    data = [[x[0], x[1],
             x[2], x[3], x[3],
             x[6], x[7], x[7],
             x[10], x[11], x[11],
             x[14], x[15], x[15],
             x[16]
             ] for x in data]
    data.extend(matches)

    t = '    Matched between exposures: %d\n' % count + \
        '    Short only: %d\n' % (len(short_data) - count) + \
        '    Long only: %d\n' % (len(data) - len(short_data)) + \
        '    Total: %d' % len(data)

    print(t)

    ax1.axhline(1, linestyle='--', color='#3d3d3d')
    ax2.axhline(0, linestyle='--', color='#3d3d3d')

    ax1.set_xlabel('Long Flux (counts/sec)')
    ax1.set_ylabel('(Short / Long) Flux Ratio')

    ax2.set_xlabel('Long Mag')
    ax2.set_ylabel('(Long-Short) Mag')

    ax1.set_ylim([0.5, 1.40])
    # ax1.set_ylim([-2000, 3000])
    ax1.set_xlim([0, 25000])

    ax2.invert_xaxis()

    ax2.set_ylim([-0.75, 0.75])

    for ax in (ax1, ax2):
        spine_lw = 4
        [ax.spines[axis].set_linewidth(spine_lw)
         for axis in ['top', 'bottom', 'left', 'right']]

        ax.legend(['B', 'V', 'R', 'H'])

    plt.tight_layout()

    directory = 'output/%s/exp_diffs/' % cluster
    if not os.path.exists(directory):
        os.makedirs(directory)

    filename = directory + '%s.png' % date
    fig.savefig(filename)

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
        bool: Returns False if point in image is not saturated, and True
              if it is saturated.

    """
    saturation = 60000

    x_floor = floor(x)
    y_floor = floor(y)

    filename = 'photometry/%s/%s/%s3.fits' % (cluster, date, fil)
    with fits.open(filename) as hdu:
        data = hdu[0].data
        for i in range(-1, 2):
            for j in range(-1, 2):
                val = data[y_floor + i, x_floor + j]
                saturated = (val >= saturation)
                if saturated:
                    break
            else:
                continue
            break

        # Return true if saturated
        return saturated


def GetRaDecs(cluster, date, data):
    """Creates a list of RA/Decs corresponding to list of matched targets.

    Args:
        cluster (str): The cluster from which the corresponding `.wcs` is
                       located.
        date (str): The date from which the corresponding `.wcs` is located.
        data (list): Data set from which the RA/Decs should be extracted

    """
    w = GetWCS(cluster, date)

    coords = []
    for target in data:
        radec = w.all_pix2world(target[0], target[1], 0)
        radec = tuple([float(i) for i in radec])
        coords.append(radec)

    filename = 'photometry/%s/%s/phot_radec.csv' % (cluster, date)
    with open(filename, 'w') as F:
        writer = csv.writer(F)
        writer.writerows(coords)


def ExtinctionCorrection(cluster, data):
    """Applies extinction correction.

    Args:
        app (Application): The GUI application object that controls processing.
        data (list): Data set on which to apply corrections.

    Returns:
        list: Data set with corrected photometry.

    """
    filename = 'photometry/%s/extinctions.dat' % cluster
    A_b, A_v, A_r = np.loadtxt(filename).tolist()

    for target in data:
        target[2] -= A_b
        target[5] -= A_v
        target[8] -= A_r

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
        filename = 'standards/%s/%s_aperture_corrections.dat' % (date, cluster)
        corrections = np.loadtxt(filename)
    except IOError:
        print("  Warning: Aperture corrections not found.")
        return data

    for target in data:
        target[2] += corrections[0]
        target[5] += corrections[1]
        target[8] += corrections[2]

    return data
