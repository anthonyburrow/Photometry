from .observations import ListDates
from .astrometry import GetAstrometryOffset
from .spectral_type import AbsMag, GetSpectralType
from .distance import GetRParams, GetDistanceOutliers
from .match import MatchTarget

import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator, FormatStrFormatter

import numpy as np
from astropy.io import fits
from astropy import wcs
from astropy.coordinates import SkyCoord
from astropy import units as u
import warnings


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
    print("\nCompiling Be list...")

    acc_count = 1
    rej_count = 1
    baseDate = ListDates(cluster)[0]

    for date in ListDates(cluster):
        filename = 'output/' + cluster + '/' + date + \
                   '/belist_scaled_' + app.phot_type + '.dat'
        data = np.loadtxt(filename)

        # Get header info
        filename = 'photometry/' + cluster + '/' + date + '/B1.fits'
        F = fits.getheader(filename)
        julian = F['JD']
        binning = F['XBINNING']
        ra_default = F['RA']
        dec_default = F['DEC']
        # Reformat ra/dec
        coo = SkyCoord(ra_default + dec_default, unit=(u.hourangle, u.deg))
        ra = coo.ra.deg
        dec = coo.dec.deg

        # Get ra/decs of CERTAIN outliers (by distance) for this date
        outliers = GetDistanceOutliers(cluster, date)

        # Read WCS file for RA and Dec transformations
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            try:
                filename = 'photometry/' + cluster + '/' + date + \
                           '/B1_wcs.fits'
                w = wcs.WCS(filename)
            except IOError:
                print("\nError: Retrieve 'B1_wcs.fits' file for " +
                      cluster + " on " + date +
                      " before calculating exact RA and Dec values. \
                      Image center values added as placeholder.")

        if date == baseDate:
            for target in data:
                try:
                    radec = w.all_pix2world(target[0], target[1], 0)
                    ra = float(radec[0])
                    dec = float(radec[1])
                except Exception:
                    ra = '--'
                    dec = '--'

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
                    if match in BeCandidates:
                        accept = True
                    else:
                        accept = False
                # If target does not match with a listed star
                else:
                    try:
                        radec = w.all_pix2world(target[0], target[1], 0)
                        ra = float(radec[0])
                        dec = float(radec[1])
                    except Exception:
                        ra = '--'
                        dec = '--'

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
        filename = 'output/' + cluster + '/' + date + \
                   '/phot_scaled_' + app.phot_type + '.dat'
        data = np.loadtxt(filename)

        # Get header info
        F = fits.getheader('photometry/' + cluster + '/' + date + '/B1.fits')
        julian = F['JD']
        binning = F['XBINNING']

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
    filename = 'output/' + cluster + '/' + date + \
               '/phot_scaled_' + app.phot_type + '.dat'
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
    F = fits.getheader('photometry/' + cluster + '/' + date + '/B1.fits')
    binning = F['XBINNING']

    for b in reversed(BData):
        xref = b[0] * binning + xOffset
        yref = b[1] * binning + yOffset
        for be in BeCandidates + rej_BeCandidates:
            if abs(xref - be[2]) <= app.cooTol and \
               abs(yref - be[3]) <= app.cooTol:
                BData.remove(b)
                break

    # Get rid of stars not in cluster (by distance)
    outliers = GetDistanceOutliers(cluster, date)

    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        try:
            w = wcs.WCS('photometry/' + cluster + '/' + date + '/B1_wcs.fits')
        except IOError:
            print("\nError: Retrieve 'B1_wcs.fits' file for " +
                  cluster + " on " + date +
                  " before calculating exact RA and Dec values. \
                  Image center values added as placeholder.")

    for target in reversed(BData):
        radec = w.all_pix2world(target[0], target[1], 0)
        ra = float(radec[0])
        dec = float(radec[1])

        if not rejected:
            # If ra/dec corresponds to distance outlier, skip target
            if [ra, dec] in outliers:
                BData.remove(target)
        else:
            # Gets distance outliers only
            if [ra, dec] not in outliers:
                BData.remove(target)

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
    filename = 'output/' + cluster + '/' + date + \
               '/beList_scaled_' + app.phot_type + '.dat'
    beData = np.loadtxt(filename)
    NumBe = len(beData)

    bData = FindBStars(cluster, app, date, BeCandidates, rej_BeCandidates)
    NumB = len(bData)

    Be_ratio = NumBe / (NumB + NumBe)

    summary_file.write("Be candidates:                " + str(NumBe) + "\n")
    summary_file.write("B stars:                      " + str(NumB) + "\n")
    summary_file.write("Be ratio:                     " +
                       "%.3f" % Be_ratio + "\n\n\n")

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
    filename = 'output/' + cluster + '/summary_' + app.phot_type + '.txt'
    summary_file = open(filename, 'w')

    t = "================================================\n" + \
        "                " + cluster + " Summary                  \n" + \
        "================================================\n\n\n"
    summary_file.write(t)

    mostB = [0, '']
    dates = ListDates(cluster)
    for date in dates:
        t = "                   " + date + "                   \n" + \
            "------------------------------------------------\n"
        summary_file.write(t)

        # Technical information
        filename = 'photometry/' + cluster + '/' + date + '/B1.fits'
        G = fits.getheader(filename)
        summary_file.write("Binning: " + str(G['XBINNING']) + '\n')

        summary_file.write("Exposures: ")
        files = ['B1', 'B3', 'V1', 'V3', 'R1', 'R3', 'H1', 'H3']
        for x in files:
            filename = 'photometry/' + cluster + '/' + date + \
                       '/' + x + '.fits'
            G = fits.getheader(filename)
            summary_file.write(x + ": " + '%d' % G['EXPTIME'])
            if x != files[len(files) - 1]:
                summary_file.write(" | ")

        summary_file.write('\n\n')

        # Threshold information
        filename = 'output/' + cluster + '/' + date + \
                   '/thresholds_' + app.phot_type + '.dat'
        thresholds = np.loadtxt(filename)
        t = "Thresholds:\n" + \
            "   Constant: " + "%.3f" % thresholds[0][1] + \
            "     Linear: " + "%.3f" % thresholds[1][0] + ", " + \
                              "%.3f" % thresholds[1][1] + "\n\n"
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

    ra = [x[5] for x in data]
    dec = [x[6] for x in data]
    julian = [x[7] for x in data]
    date = [x[8] for x in data]

    count = [x[9] for x in data]
    for i in range(1, max(count) + 1):
        if i not in count:
            count = [x - 1 for x in count if x > i]

    bmag = [x[10] for x in data]
    berr = [x[11] for x in data]
    vmag = [x[12] for x in data]
    verr = [x[13] for x in data]
    rmag = [x[14] for x in data]
    rerr = [x[15] for x in data]
    hmag = [x[16] for x in data]
    herr = [x[17] for x in data]

    distance_mean, distance_std = GetRParams(cluster)

    absVmag = [AbsMag(x, distance_mean) for x in vmag]   # 19
    spectralTypes = [GetSpectralType(x, distance_mean) for x in vmag]   # 20
    for i in range(len(data)):
        data[i].extend([absVmag[i], spectralTypes[i]])

    # Open summary files
    if belist == BeCandidates:
        filename = 'output/' + cluster + \
                   '/BeList_' + app.phot_type + '.txt'
    else:
        filename = 'output/' + cluster + \
                   '/rejected_BeList_' + app.phot_type + '.txt'
    list_file = open(filename, 'w')

    for i in range(len(data)):
        # Get numerical excess
        try:
            filename = 'output/' + cluster + '/' + date[i] + \
                       '/thresholds_' + app.phot_type + '.dat'
            thresholds = np.loadtxt(filename)
            if app.threshold_type == 'Constant':
                slope = thresholds[0][0]
                intercept = thresholds[0][1]
            elif app.threshold_type == 'Linear':
                slope = thresholds[1][0]
                intercept = thresholds[1][1]
        except Exception:
            print("Error: Threshold calculations are required.")

        r_h = rmag[i] - hmag[i]
        b_v = bmag[i] - vmag[i]
        r_h0 = slope * b_v + intercept
        excess = r_h - r_h0

        # Write to file
        t = cluster + '-WBBe' + str(count[i]) + '\t' + \
            '%.3f' % absVmag[i] + '\t' + \
            spectralTypes[i] + '\t' + \
            '%.10f' % ra[i] + '\t' + \
            '%.10f' % dec[i] + '\t' + \
            '%.10f' % julian[i] + '\t' + \
            '%.3f' % bmag[i] + '\t' + \
            '%.3f' % berr[i] + '\t' + \
            '%.3f' % vmag[i] + '\t' + \
            '%.3f' % verr[i] + '\t' + \
            '%.3f' % rmag[i] + '\t' + \
            '%.3f' % rerr[i] + '\t' + \
            '%.3f' % hmag[i] + '\t' + \
            '%.3f' % herr[i] + '\t' + \
            '%.3f' % excess + '\n'

        list_file.write(t)

    list_file.close()

    if belist != BeCandidates:
        return

    be_spectralTypes = []
    for i in range(max(count)):
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

    filename = 'output/' + cluster + \
               '/spectral_types_' + app.phot_type + '.png'
    fig.savefig(filename)
    plt.clf()

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
            "    Be: " + str(len(be_type_unknown)) + '\n' + \
            "    B: " + str(len(b_type_unknown)) + '\n' + \
            "    Ratio: " + '%.3f' % get_ratio(be_type_unknown, b_type_unknown) + \
            '\n\n' + \
            "Type O:\n" + \
            "    Be: " + str(len(be_type_O)) + '\n' + \
            "    B: " + str(len(b_type_O)) + '\n' + \
            "    Ratio: " + '%.3f' % get_ratio(be_type_O, b_type_O) + \
            '\n\n' + \
            "Type B0-B3:\n" + \
            "    Be: " + str(len(be_type_B0_B3)) + '\n' + \
            "    B: " + str(len(b_type_B0_B3)) + '\n' + \
            "    Ratio: " + '%.3f' % get_ratio(be_type_B0_B3, b_type_B0_B3) + \
            '\n\n' + \
            "Type B4-B5:\n" + \
            "    Be: " + str(len(be_type_B4_B5)) + '\n' + \
            "    B: " + str(len(b_type_B4_B5)) + '\n' + \
            "    Ratio: " + '%.3f' % get_ratio(be_type_B4_B5, b_type_B4_B5) + \
            '\n\n' + \
            "Type B6-B9:\n" + \
            "    Be: " + str(len(be_type_B6_B9)) + '\n' + \
            "    B: " + str(len(b_type_B6_B9)) + '\n' + \
            "    Ratio: " + '%.3f' % get_ratio(be_type_B6_B9, b_type_B6_B9) + \
            '\n\n' + \
            "Type A:\n" + \
            "    Be: " + str(len(be_type_A)) + '\n' + \
            "    B: " + str(len(b_type_A)) + '\n' + \
            "    Ratio: " + '%.3f' % get_ratio(be_type_A, b_type_A) + \
            '\n\n'

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
    filename = 'output/' + cluster + '/' + date + \
               '/belist_scaled_' + app.phot_type + '.dat'
    Be_all = np.loadtxt(filename).tolist()

    # Be candidates in cluster
    Be_in = []
    for target in Be_all:
        for candidate in BeCandidates:
            if target[0] == candidate[0] and target[1] == candidate[1] and \
               date == candidate[8]:
                Be_in.append(target)

    # Be candidates outside cluster
    Be_out = []
    for target in Be_all:
        for candidate in rej_BeCandidates:
            if target[0] == candidate[0] and target[1] == candidate[1] and \
               date == candidate[8]:
                Be_out.append(target)

    # B-type stars in cluster
    B_in = FindBStars(cluster, app, date, BeCandidates, rej_BeCandidates)

    # B-type stars outside cluster
    B_out = FindBStars(cluster, app, date, BeCandidates, rej_BeCandidates,
                       rejected=True)

    # Plot 2CD
    filename = 'output/' + cluster + '/' + date + \
               '/phot_scaled_' + app.phot_type + '.dat'
    data = np.loadtxt(filename)

    B_V = data[:, 2] - data[:, 4]
    B_Verr = np.sqrt(data[:, 3]**2 + data[:, 5]**2)
    R_H = data[:, 6] - data[:, 8]
    R_Herr = np.sqrt(data[:, 7]**2 + data[:, 9]**2)

    plt.style.use('researchpaper')
    fig, ax = plt.subplots()

    ax.plot(B_V, R_H, 'o', color='#3f3f3f', markersize=12)
    ax.set_xlabel('B-V', fontsize=36)
    ax.set_ylabel('R-Halpha', fontsize=36)

    apCorr = np.loadtxt('standards/' + date + '/' + cluster +
                        '_aperture_corrections.dat')

    ax.set_xlim([app.B_VMin - 0.1, app.B_VMax + 2.5])
    ax.set_ylim([-6.5 + apCorr[2], -4 + apCorr[2]])

    ax.xaxis.set_major_locator(MultipleLocator(1))
    ax.xaxis.set_major_formatter(FormatStrFormatter('%d'))
    ax.xaxis.set_minor_locator(MultipleLocator(0.25))

    ax.yaxis.set_major_locator(MultipleLocator(1))
    ax.yaxis.set_major_formatter(FormatStrFormatter('%d'))
    ax.yaxis.set_minor_locator(MultipleLocator(0.25))

    spine_lw = 4
    [ax.spines[axis].set_linewidth(spine_lw)
     for axis in ['top', 'bottom', 'left', 'right']]

    ax.errorbar(B_V, R_H, xerr=B_Verr, yerr=R_Herr, fmt='none',
                ecolor='#50a0e5', elinewidth=7)

    # Overplot Be and B stars
    ax.plot([x[2] - x[4] for x in B_in], [x[6] - x[8] for x in B_in], 'o',
            color='#ff5151', markersize=8, markeredgewidth=5,
            label='B in cluster')
    ax.plot([x[2] - x[4] for x in B_out], [x[6] - x[8] for x in B_out], 'o',
            color='#6ba3ff', markersize=8, markeredgewidth=5,
            label='B outside cluster')
    ax.plot([x[2] - x[4] for x in Be_in], [x[6] - x[8] for x in Be_in], 'x',
            color='#e52424', markersize=15, markeredgewidth=5,
            label='Be in cluster')
    ax.plot([x[2] - x[4] for x in Be_out], [x[6] - x[8] for x in Be_out], 'x',
            color='#2f72e0', markersize=15, markeredgewidth=5,
            label='Be outside cluster')

    # Plot threshold line
    filename = 'output/' + cluster + '/' + date + \
               '/thresholds_' + app.phot_type + '.dat'
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
        ax.plot(linex, liney, '--', color='#ff5151', label='Be Threshold',
                linewidth=6)
    except IOError:
        print("\nNote: Thresholds have not been calculated or written \
               to file yet and will not be displayed.")

    ax.legend()

    # Output
    filename = 'output/' + cluster + '/' + date + \
               '/plots/2CD_' + app.phot_type + '_detailed.png'
    fig.savefig(filename)

    plt.clf()

    # Plot CMD
    plt.style.use('researchpaper')
    fig, ax = plt.subplots()

    ax.plot(B_V, data[:, 4], 'o', color='#3f3f3f', markersize=12)
    ax.set_xlabel('B-V')
    ax.set_ylabel('V')

    ax.set_xlim([app.B_VMin - 0.1, app.B_VMax + 2.5])
    ax.set_ylim([18.5 - app.A_v, 8.5 - app.A_v])

    ax.xaxis.set_major_locator(MultipleLocator(1))
    ax.xaxis.set_major_formatter(FormatStrFormatter('%d'))
    ax.xaxis.set_minor_locator(MultipleLocator(0.25))

    ax.yaxis.set_major_locator(MultipleLocator(2))
    ax.yaxis.set_major_formatter(FormatStrFormatter('%d'))
    ax.yaxis.set_minor_locator(MultipleLocator(0.5))

    spine_lw = 4
    [ax.spines[axis].set_linewidth(spine_lw)
     for axis in ['top', 'bottom', 'left', 'right']]

    ax.errorbar(B_V, data[:, 4], xerr=B_Verr, yerr=data[:, 5], fmt='none',
                ecolor='#50a0e5', elinewidth=7)

    # Overplot Be and B stars
    ax.plot([x[2] - x[4] for x in B_in], [x[4] for x in B_in], 'D',
            color='#ff5151', markersize=8, markeredgewidth=5,
            label='B in cluster')
    ax.plot([x[2] - x[4] for x in B_out], [x[4] for x in B_out], 'D',
            color='#6ba3ff', markersize=8, markeredgewidth=5,
            label='B outside cluster')
    ax.plot([x[2] - x[4] for x in Be_in], [x[4] for x in Be_in], 'x',
            color='#e52424', markersize=15, markeredgewidth=5,
            label='Be in cluster')
    ax.plot([x[2] - x[4] for x in Be_out], [x[4] for x in Be_out], 'x',
            color='#2f72e0', markersize=15, markeredgewidth=5,
            label='Be outside cluster')

    ax.legend()

    # Output
    filename = 'output/' + cluster + '/' + date + \
               '/plots/CMD_' + app.phot_type + '_detailed.png'
    fig.savefig(filename)

    plt.clf()
