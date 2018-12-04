from .observations import ListDates
from .astrometry import GetAstrometryOffset
from .spectral_type import AbsMag, GetSpectralType
from .distance import GetRParams, GetDistanceOutliers
from .match import MatchTarget
from .read_files import Binning, GetWCS, GetFITSValues

import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator, FormatStrFormatter

import numpy as np
from astropy.io import fits
# from astropy.coordinates import SkyCoord
# from astropy import units as u


def ProcessAnalysis(cluster, app):
    """Controls full analysis process.

    This program compiles a full list of Be candidates from the calculated Be
    candidates of individual nights, finds corresponding targets from other
    dates, then summarizes the results in output.

    Args:
        cluster (str): Cluster for which analysis will be conducted.
        app (Application): The GUI application object that controls processing.

    """
    # Compile list of Be stars
    BeCandidates = []
    rej_BeCandidates = []
    CompileBeLists(cluster, app, BeCandidates, rej_BeCandidates)
    FindCorrespondingTargets(cluster, app, BeCandidates)
    FindCorrespondingTargets(cluster, app, rej_BeCandidates)

    # Write to files
    mostB_date = NightSummary(cluster, app, BeCandidates, rej_BeCandidates)
    BeSummary(cluster, app, BeCandidates, mostB_date,
              BeCandidates, rej_BeCandidates)
    BeSummary(cluster, app, rej_BeCandidates, mostB_date,
              BeCandidates, rej_BeCandidates)

    # Make plots showing Be/B and accepted/rejected separation
    for date in ListDates(cluster):
        BeCandidatePlots(cluster, app, date, BeCandidates, rej_BeCandidates)


def CompileBeLists(cluster, app, BeCandidates, rej_BeCandidates):
    """Creates data sets for all Be candidates from each date of observation.

    This reads from output from each night to compile each indivual Be candidate
    and uniquely identify them.  Targets are then put into an 'accepted' or
    'rejected' list based on calculated cluster membership.

    Args:
        cluster (str): Cluster for which Be candidates will be compiled.
        app (Application): The GUI application object that controls processing.
        BeCandidates (list): List containing accepted Be candidates.
        rej_BeCandidates (list): List containing rejected Be candidates.

    """
    acc_count = 1
    rej_count = 1
    baseDate = ListDates(cluster)[0]

    for date in ListDates(cluster):
        print("  Retrieving Be data from %s..." % date)

        filename = 'output/' + cluster + '/' + date + '/belist_scaled.dat'
        data = np.loadtxt(filename, ndmin=2)

        # Get header info
        fits = GetFITSValues(cluster, date, ['JD', 'XBINNING', 'RA', 'DEC'])
        julian = fits['JD']
        binning = fits['XBINNING']

        # Set default ra/decs
        # ra_default = fits['RA']
        # dec_default = fits['DEC']
        # Reformat ra/dec
        # coo = SkyCoord(ra_default + dec_default, unit=(u.hourangle, u.deg))
        # ra = coo.ra.deg
        # dec = coo.dec.deg

        # Get ra/decs of CERTAIN outliers (by distance) for this date
        outliers = GetDistanceOutliers(cluster, date)
        w = GetWCS(cluster, date)

        if date == baseDate:
            for target in data:
                radec = w.all_pix2world(target[0], target[1], 0)
                ra, dec = tuple([float(i) for i in radec])

                if [ra, dec] not in outliers:
                    count = acc_count
                    acc_count += 1
                else:
                    count = rej_count
                    rej_count += 1

                be = [
                    target[0],                       # 0 - x on image
                    target[1],                       # 1 - y on image
                    target[0] * binning,             # 2 - x ref
                    target[1] * binning,             # 3 - y ref
                    'present',                       # 4 - transient status
                    ra,                              # 5 - ra
                    dec,                             # 6 - dec
                    julian,                          # 7 - julian date
                    date,                            # 8 - date
                    count,                           # 9 - count
                    # 10-17 - photometry
                    target[2], target[3], target[4], target[5],
                    target[6], target[7], target[8], target[9],
                    True                             # 18 - be star check
                ]

                if [ra, dec] not in outliers:
                    BeCandidates.append(be)
                else:
                    rej_BeCandidates.append(be)
        else:
            # Get x- and y-offsets
            xOffset, yOffset = GetAstrometryOffset(cluster, date, baseDate)

            # Compare with other dates, see if target is already found
            for target in data:
                xRef = binning * target[0] + xOffset
                yRef = binning * target[1] + yOffset
                # If target refers to an already listed star
                match = MatchTarget(app.cooTol, [xRef, yRef],
                                    BeCandidates + rej_BeCandidates, [2, 3])
                if match is not None:
                    ra = match[5]
                    dec = match[6]
                    count = match[9]
                    accept = match in BeCandidates

                # If target does not match with a listed star
                else:
                    radec = w.all_pix2world(target[0], target[1], 0)
                    ra, dec = tuple([float(i) for i in radec])

                    if [ra, dec] not in outliers:
                        count = acc_count
                        acc_count += 1
                        accept = True
                    else:
                        count = rej_count
                        rej_count += 1
                        accept = False

                be = [
                    target[0],                       # 0 - x on image
                    target[1],                       # 1 - y on image
                    xRef,                            # 2 - x ref
                    yRef,                            # 3 - y ref
                    'present',                       # 4 - transient status
                    ra,                              # 5 - ra
                    dec,                             # 6 - dec
                    julian,                          # 7 - julian date
                    date,                            # 8 - date
                    count,                           # 9 - count
                    # 10-17 - photometry
                    target[2], target[3], target[4], target[5],
                    target[6], target[7], target[8], target[9],
                    True                             # 18 - be star check
                ]

                if accept:
                    BeCandidates.append(be)
                else:
                    rej_BeCandidates.append(be)


def FindCorrespondingTargets(cluster, app, belist):
    """Finds targets from other nights that correspond to Be candidates.

    Uses matching to determine which targets from each night of observation
    identify with each star in given Be candidate list.

    Args:
        cluster (str): Cluster for which Be candidates will be compiled.
        app (Application): The GUI application object that controls processing.
        belist (list): List (accepted or rejected) to check for corresponding
                       targets.

    """
    dates = ListDates(cluster)
    baseDate = dates[0]

    targets_to_lookup = list(belist)

    for date in dates:
        filename = 'output/' + cluster + '/' + date + '/phot_scaled.dat'
        data = np.loadtxt(filename)

        # Get header info
        fits = GetFITSValues(cluster, date, ['JD', 'XBINNING'])
        julian = fits['JD']
        binning = fits['XBINNING']

        # Get x- and y-offsets
        if date != baseDate:
            xOffset, yOffset = GetAstrometryOffset(cluster, date, baseDate)
        else:
            xOffset = 0
            yOffset = 0

        _targets_to_lookup = list(targets_to_lookup)

        for i in range(1, max([x[9] for x in targets_to_lookup]) + 1):
            candidate = [x for x in _targets_to_lookup if x[9] == i]
            if any(x[8] == date for x in candidate):
                _targets_to_lookup = [x for x in _targets_to_lookup
                                      if x not in candidate]

        for target in data:
            xRef = binning * target[0] + xOffset
            yRef = binning * target[1] + yOffset
            match = MatchTarget(app.cooTol, [xRef, yRef],
                                _targets_to_lookup, [2, 3])
            if match is None:
                continue

            be = [
                target[0],                         # 0 - x on image
                target[1],                         # 1 - y on image
                xRef,                              # 2 - x ref
                yRef,                              # 3 - y ref
                'absent',                          # 4 - transient status
                match[5],                          # 5 - ra
                match[6],                          # 6 - dec
                julian,                            # 7 - julian date
                date,                              # 8 - date
                match[9],                          # 9 - count
                # 10-17 - photometry
                target[2], target[3], target[4], target[5],
                target[6], target[7], target[8], target[9],
                False                            # 18 - be star check
            ]
            belist.append(be)
            _targets_to_lookup.remove(match)


def FindBStars(cluster, app, date, BeCandidates, rej_BeCandidates,
               rejected=False):
    """Finds B-type stars without H-alpha excess for given date.

    Using photometric bounds and calculated Be candidates as well as cluster
    membership calculations, B-type targets are extracted.

    Args:
        cluster (str): Cluster for which B stars will be found.
        app (Application): The GUI application object that controls processing.
        date (str): Date for which B stars will be found.
        BeCandidates (list): List of accepted Be candidates.
        rej_BeCandidates (list): List of rejected Be candidates.
        rejected (bool): Determines if accepted or rejected B-type stars are
                         desired. False if accepted is desired, and vice versa.

    Returns:
        list: List of B-type targets

    """
    if not rejected:
        filename = 'output/' + cluster + '/' + date + '/phot_scaled_accepted.dat'
    else:
        filename = 'output/' + cluster + '/' + date + '/phot_scaled_rejected.dat'

    data = np.loadtxt(filename).tolist()
    BData = []

    # Find all B- and Be-type stars
    for target in data:
        b_v = target[2] - target[4]
        if app.B_VMin < b_v < app.B_VMax and target[4] <= 13.51:
            BData.append(target)

    # Pick out observed Be stars to result in ONLY B-type stars
    baseDate = ListDates(cluster)[0]
    xOffset, yOffset = GetAstrometryOffset(cluster, date, baseDate)
    binning = Binning(cluster, date)

    for b in reversed(BData):
        xref = b[0] * binning + xOffset
        yref = b[1] * binning + yOffset
        for be in BeCandidates + rej_BeCandidates:
            if abs(xref - be[2]) <= app.cooTol and \
               abs(yref - be[3]) <= app.cooTol:
                BData.remove(b)
                break

    return BData


def BeValues(cluster, app, date, summary_file, mostB, BeCandidates,
             rej_BeCandidates):
    """Writes numbers and ratios into a summary file for a single date.

    Args:
        cluster (str): Cluster for which information is written.
        app (Application): The GUI application object that controls processing.
        date (str): Date for which information is written.
        summary_file (file): File in which information is stored.
        mostB (list): Previously determined max number of B-type stars and
                      corresponding date.
        BeCandidates (list): List of accepted Be candidates.
        rej_BeCandidates (list): List of rejected Be candidates.

    Returns:
        list: List containing the number of B-type targets cumulatively for the
              night with the most B-type stars along with the date itself

    """
    filename = 'output/' + cluster + '/' + date + '/beList_scaled.dat'
    beData = np.loadtxt(filename)
    NumBe = len(beData)

    bData = FindBStars(cluster, app, date, BeCandidates, rej_BeCandidates)
    NumB = len(bData)

    try:
        Be_ratio = NumBe / (NumB + NumBe)
    except ZeroDivisionError:
        Be_ratio = 0

    summary_file.write("Be candidates:              %d\n" % NumBe)
    summary_file.write("B stars:                    %d\n" % NumB)
    summary_file.write("Be ratio:                   %.3f\n\n\n" % Be_ratio)

    mostB_new = mostB
    if NumB > mostB[0]:
        mostB_new = [NumB, date]

    return mostB_new


def NightSummary(cluster, app, BeCandidates, rej_BeCandidates):
    """Writes a summary of each night to a file.

    Args:
        cluster (str): Cluster for which information is written.
        app (Application): The GUI application object that controls processing.
        BeCandidates (list): List of accepted Be candidates.
        rej_BeCandidates (list): List of rejected Be candidates.

    Returns:
        str: Returns the date with the most B-type targets observed.

    """
    filename = 'output/' + cluster + '/summary.txt'
    summary_file = open(filename, 'w')

    t = '================================================\n' + \
        '                %s Summary                  \n' % cluster + \
        '================================================\n\n\n'
    summary_file.write(t)

    mostB = [0, '']
    dates = ListDates(cluster)
    for date in dates:
        t = '                   %s                   \n' % date + \
            '------------------------------------------------\n'
        summary_file.write(t)

        # Technical information
        filename = 'photometry/' + cluster + '/' + date + '/B1.fits'
        binning = Binning(cluster, date)
        summary_file.write('Binning: %d\n' % binning)

        summary_file.write('Exposures: ')
        files = ['B1', 'B3', 'V1', 'V3', 'R1', 'R3', 'H1', 'H3']
        for x in files:
            filename = 'photometry/' + cluster + '/' + date + \
                       '/' + x + '.fits'
            G = fits.getheader(filename)
            summary_file.write('%s: %d' % (x, G['EXPTIME']))
            if x != files[len(files) - 1]:
                summary_file.write(' | ')

        summary_file.write('\n\n')

        # Threshold information
        filename = 'output/' + cluster + '/' + date + '/thresholds.dat'
        thresholds = np.loadtxt(filename)
        t = "Thresholds:\n" + \
            "   Constant: %.3f     Linear: %.3f, %.3f\n\n" % \
            (thresholds[0][1], thresholds[1][0], thresholds[1][1])

        summary_file.write(t)

        # Be candidate information
        mostB = BeValues(cluster, app, date, summary_file, mostB,
                         BeCandidates, rej_BeCandidates)

    summary_file.close()

    return mostB[1]


def BeSummary(cluster, app, belist, mostB_date, BeCandidates,
              rej_BeCandidates):
    """Creates finalized output for Be candidate lists.

    This finalized summary includes the determination of each target's spectral
    type and its excess R-H value from the calculated threshold for each night.
    Also analyzes spectral type distribution through histograms and ratios.

    Args:
        cluster (str): Cluster for which information is written.
        app (Application): The GUI application object that controls processing.
        belist (list): Specific list of Be candidates to analyze.
        mostB_date (str): Date with the most observed B-type targets.
        BeCandidates (list): List of accepted Be candidates.
        rej_BeCandidates (list): List of rejected Be candidates.

    """
    if not belist:
        return

    # Sort by identifier (count) then date
    data = sorted(belist, key=lambda x: (x[9], x[8]))

    count = [x[9] for x in data]
    for i in range(1, max(count) + 1):
        if i not in count:
            count = [x - 1 for x in count if x > i]
    for i in range(len(data)):
        data[i][9] = count[i]

    # Give each target an abs. magnitude & spectral type
    distance_mean, distance_std = GetRParams(cluster)

    for target in data:
        target.extend([AbsMag(target[12], distance_mean),
                       GetSpectralType(target[12], distance_mean)])

    # Open summary files
    if belist == BeCandidates:
        filename = 'output/' + cluster + '/BeList.txt'
    else:
        filename = 'output/' + cluster + '/rejected_BeList.txt'
    list_file = open(filename, 'w')

    for target in data:
        # Get numerical excess
        filename = 'output/' + cluster + '/' + target[8] + \
                   '/thresholds.dat'
        thresholds = np.loadtxt(filename)
        if app.threshold_type == 'Constant':
            slope = thresholds[0][0]
            intercept = thresholds[0][1]
        elif app.threshold_type == 'Linear':
            slope = thresholds[1][0]
            intercept = thresholds[1][1]

        r_h = target[14] - target[16]
        b_v = target[10] - target[12]
        r_h0 = slope * b_v + intercept
        excess = r_h - r_h0

        # Write to file
        t = '%s-WBBe%d' % (cluster, target[9]) + '\t' + \
            '%.3f' % target[19] + '\t' + \
            target[20] + '\t' + \
            '%.10f' % target[5] + '\t' + \
            '%.10f' % target[6] + '\t' + \
            '%.10f' % target[7] + '\t' + \
            '%.3f' % target[10] + '\t' + \
            '%.3f' % target[11] + '\t' + \
            '%.3f' % target[12] + '\t' + \
            '%.3f' % target[13] + '\t' + \
            '%.3f' % target[14] + '\t' + \
            '%.3f' % target[15] + '\t' + \
            '%.3f' % target[16] + '\t' + \
            '%.3f' % target[17] + '\t' + \
            '%.3f' % excess + '\n'

        list_file.write(t)

    list_file.close()

    if belist != BeCandidates:
        return

    be_spectralTypes = []
    for i in range(min(count), max(count)):
        m = [x[12] for x in data if x[9] == i]
        s = GetSpectralType(np.mean(m), distance_mean)
        be_spectralTypes.append(s)

    be_type_unknown = [x for x in be_spectralTypes if x == '--']
    be_type_O = [x for x in be_spectralTypes if x[0] == 'O']
    be_type_B0_B3 = [x for x in be_spectralTypes
                     if x == 'B0' or x == 'B1' or x == 'B2' or x == 'B3']
    be_type_B4_B5 = [x for x in be_spectralTypes
                     if x == 'B4' or x == 'B5']
    be_type_B6_B9 = [x for x in be_spectralTypes
                     if x == 'B6' or x == 'B7' or x == 'B8' or x == 'B9']
    be_type_A = [x for x in be_spectralTypes if x[0] == 'A']

    frequencies = [len(be_type_unknown), len(be_type_O), len(be_type_B0_B3),
                   len(be_type_B4_B5), len(be_type_B6_B9), len(be_type_A)]
    names = ['<O6', 'O6-O8', 'B0-B3', 'B4-B5', 'B6-B9', 'A']

    # Plot spectral type histogram
    plt.style.use('researchpaper')
    fig, ax = plt.subplots()

    x_coordinates = np.arange(len(frequencies))
    ax.bar(x_coordinates, frequencies, align='center', color='#3f3f3f')
    ax.set_xlabel('Spectral Type')
    ax.set_ylabel('Frequency')

    ax.xaxis.set_major_locator(plt.FixedLocator(x_coordinates))
    ax.xaxis.set_major_formatter(plt.FixedFormatter(names))

    spine_lw = 4
    [ax.spines[axis].set_linewidth(spine_lw)
     for axis in ['top', 'bottom', 'left', 'right']]

    filename = 'output/' + cluster + '/spectral_types.png'
    fig.savefig(filename)

    plt.close('all')

    # Be ratios by spectral type
    BData = FindBStars(cluster, app, mostB_date, BeCandidates,
                       rej_BeCandidates)
    b_spectralTypes = []
    for target in BData:
        b_spectralTypes.append(GetSpectralType(target[4], distance_mean))

    b_type_unknown = [x for x in b_spectralTypes if x[0] == '--']
    b_type_O = [x for x in b_spectralTypes if x[0] == 'O']
    b_type_B0_B3 = [x for x in b_spectralTypes
                    if x == 'B0' or x == 'B1' or x == 'B2' or x == 'B3']
    b_type_B4_B5 = [x for x in b_spectralTypes if x == 'B4' or x == 'B5']
    b_type_B6_B9 = [x for x in b_spectralTypes
                    if x == 'B6' or x == 'B7' or x == 'B8' or x == 'B9']
    b_type_A = [x for x in b_spectralTypes if x[0] == 'A']

    def get_ratio(be, b):
        try:
            ratio = len(be) / (len(be) + len(b))
        except ZeroDivisionError:
            ratio = 0
        return ratio

    filename = 'output/' + cluster + '/ratios_by_spec_type.txt'
    with open(filename, 'w') as F:
        t = "Unknown:\n" + \
            '    Be: %d\n' % len(be_type_unknown) + \
            '    B: %d\n' % len(b_type_unknown) + \
            '    Ratio: %.3f\n\n' % get_ratio(be_type_unknown, b_type_unknown) + \
            'Type O:\n' + \
            '    Be: %d\n' % len(be_type_O) + \
            '    B: %d\n' % len(b_type_O) + \
            '    Ratio: %.3f\n\n' % get_ratio(be_type_O, b_type_O) + \
            'Type B0-B3:\n' + \
            '    Be: %d\n' % len(be_type_B0_B3) + \
            '    B: %d\n' % len(b_type_B0_B3) + \
            '    Ratio: %.3f\n\n' % get_ratio(be_type_B0_B3, b_type_B0_B3) + \
            'Type B4-B5:\n' + \
            '    Be: %d\n' % len(be_type_B4_B5) + \
            '    B: %d\n' % len(b_type_B4_B5) + \
            '    Ratio: %.3f\n\n' % get_ratio(be_type_B4_B5, b_type_B4_B5) + \
            'Type B6-B9:\n' + \
            '    Be: %d\n' % len(be_type_B6_B9) + \
            '    B: %d\n' % len(b_type_B6_B9) + \
            '    Ratio: %.3f\n\n' % get_ratio(be_type_B6_B9, b_type_B6_B9) + \
            'Type A:\n' + \
            '    Be: %d\n' % len(be_type_A) + \
            '    B: %d\n' % len(b_type_A) + \
            '    Ratio: %.3f\n\n' % get_ratio(be_type_A, b_type_A)

        F.write(t)


def BeCandidatePlots(cluster, app, date, BeCandidates, rej_BeCandidates):
    """Creates detailed CMDs and 2CDs for a specified date.

    Detailed plots created are CMDs and 2CDs that specify cluster membership
    of all observed Be candidates and B-type stars.

    Args:
        cluster (str): Cluster for which information is written.
        app (Application): The GUI application object that controls processing.
        date (str): Date for which information is written.
        BeCandidates (list): List of accepted Be candidates.
        rej_BeCandidates (list): List of rejected Be candidates.

    """
    filename = 'output/' + cluster + '/' + date + '/phot_scaled.dat'
    phot = np.loadtxt(filename, ndmin=2).tolist()

    # Be candidates in cluster
    Be_in = []
    for candidate in (x for x in BeCandidates if x[8] == date):
        for target in phot:
            if target[0] == candidate[0] and target[1] == candidate[1]:
                Be_in.append(target)
                phot.remove(target)
                break

    # Be candidates outside cluster
    Be_out = []
    for candidate in (x for x in rej_BeCandidates if x[8] == date):
        for target in phot:
            if target[0] == candidate[0] and target[1] == candidate[1]:
                Be_out.append(target)
                phot.remove(target)
                break

    # B-type stars in cluster
    B_in = FindBStars(cluster, app, date, BeCandidates, rej_BeCandidates)

    # B-type stars outside cluster
    B_out = FindBStars(cluster, app, date, BeCandidates, rej_BeCandidates,
                       rejected=True)

    # Plot data
    filename = 'output/' + cluster + '/' + date + '/phot_scaled.dat'
    data = np.loadtxt(filename)

    B_V = data[:, 2] - data[:, 4]
    B_Verr = np.sqrt(data[:, 3]**2 + data[:, 5]**2)
    R_H = data[:, 6] - data[:, 8]
    R_Herr = np.sqrt(data[:, 7]**2 + data[:, 9]**2)

    def Plot(plot_type):
        plt.style.use('researchpaper')
        fig, ax = plt.subplots()

        apCorr = np.loadtxt('standards/' + date + '/' + cluster +
                            '_aperture_corrections.dat')

        if plot_type == 'cmd':
            # Plotting
            ax.plot(B_V, data[:, 4], 'o', color='#3d3d3d', markersize=10)
            ax.errorbar(B_V, data[:, 4], xerr=B_Verr, yerr=data[:, 5],
                        fmt='none', ecolor='#8c8c8c', elinewidth=7)

            ax.plot([x[2] - x[4] for x in B_out], [x[4] for x in B_out], 'o',
                    color='#6ba3ff', markersize=12, markeredgewidth=2,
                    label='B outside cluster')
            ax.plot([x[2] - x[4] for x in B_in], [x[4] for x in B_in], 'o',
                    color='#ff5151', markersize=12, markeredgewidth=2,
                    label='B in cluster')
            ax.plot([x[2] - x[4] for x in Be_out], [x[4] for x in Be_out], 'D',
                    color='#6ba3ff', markersize=12, markeredgewidth=2,
                    markeredgecolor='#2d2d2d', label='Be outside cluster')
            ax.plot([x[2] - x[4] for x in Be_in], [x[4] for x in Be_in], 'D',
                    color='#ff5151', markersize=12, markeredgewidth=2,
                    markeredgecolor='#2d2d2d', label='Be in cluster')

            # Settings
            ax.set_ylabel('V')
            ax.set_ylim([18.5 - app.A_v, 8.5 - app.A_v])

            ax.yaxis.set_major_locator(MultipleLocator(2))
            ax.yaxis.set_major_formatter(FormatStrFormatter('%d'))
            ax.yaxis.set_minor_locator(MultipleLocator(0.5))

            output = 'CMD_detailed.png'

        elif plot_type == '2cd':
            # Plotting
            ax.plot(B_V, R_H, 'o', color='#3d3d3d', markersize=10)
            ax.errorbar(B_V, R_H, xerr=B_Verr, yerr=R_Herr, fmt='none',
                        ecolor='#8c8c8c', elinewidth=7)

            ax.plot([x[2] - x[4] for x in B_out], [x[6] - x[8] for x in B_out],
                    'o', color='#6ba3ff', markersize=12, markeredgewidth=2,
                    label='B outside cluster')
            ax.plot([x[2] - x[4] for x in B_in], [x[6] - x[8] for x in B_in],
                    'o', color='#ff5151', markersize=12, markeredgewidth=2,
                    label='B in cluster')
            ax.plot([x[2] - x[4] for x in Be_out], [x[6] - x[8] for x in Be_out],
                    'D', color='#6ba3ff', markersize=12, markeredgewidth=2,
                    markeredgecolor='#2d2d2d', label='Be outside cluster')
            ax.plot([x[2] - x[4] for x in Be_in], [x[6] - x[8] for x in Be_in],
                    'D', color='#ff5151', markersize=12, markeredgewidth=2,
                    markeredgecolor='#2d2d2d', label='Be in cluster')

            # Settings
            ax.set_ylabel('R-Halpha')
            ax.set_ylim([-6.5 + apCorr[2], -4 + apCorr[2]])

            ax.yaxis.set_major_locator(MultipleLocator(1))
            ax.yaxis.set_major_formatter(FormatStrFormatter('%d'))
            ax.yaxis.set_minor_locator(MultipleLocator(0.25))

            filename = 'output/' + cluster + '/' + date + '/thresholds.dat'
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

            output = '2CD_detailed.png'

        # Shared settings
        ax.set_xlabel('B-V')
        ax.set_xlim([app.B_VMin - 0.1, app.B_VMax + 2.5])

        ax.xaxis.set_major_locator(MultipleLocator(1))
        ax.xaxis.set_major_formatter(FormatStrFormatter('%d'))
        ax.xaxis.set_minor_locator(MultipleLocator(0.25))

        spine_lw = 4
        [ax.spines[axis].set_linewidth(spine_lw)
         for axis in ['top', 'bottom', 'left', 'right']]

        ax.legend()

        # Output
        filename = 'output/' + cluster + '/' + date + \
                   '/plots/' + output
        fig.savefig(filename)

        plt.close('all')

    Plot('2cd')
    Plot('cmd')
