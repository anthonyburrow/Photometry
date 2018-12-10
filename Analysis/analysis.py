from .observations import ListDates
from .astrometry import GetAstrometryOffset
from .spectral_type import AbsMag, SpectralType
from .distance import GetRParams, GetDistanceOutliers
from .read_files import Binning, GetWCS, GetFITSValues

import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator, FormatStrFormatter

import numpy as np
from astropy.io import fits


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


class BeCandidate:

    def __init__(self):
        self.x = None
        self.y = None
        self.xref = None
        self.yref = None
        self.ra = None
        self.dec = None
        self.julian = None
        self.date = None
        self.count = None

        self.phot = {}

        self.absmag = None
        self.spectraltype = None

        self.observed_be = True

    def GetAbsMag(self, distance):
        if self.absmag is None:
            self.absmag = AbsMag(self.phot['Vmag'], distance)

        return self.absmag

    def GetSpectralType(self, distance):
        if self.spectraltype is None:
            self.spectraltype = SpectralType(self.phot['Vmag'], distance)

        return self.spectraltype


def MatchTarget(tol, coords, data):
    potentialMatches = []
    for target in data:
        if abs(coords[0] - target.xref) <= tol and \
           abs(coords[1] - target.yref) <= tol:
            potentialMatches.append(target)

    # Match with the closest target
    if potentialMatches:
        d = [np.sqrt((coords[0] - x.xref)**2 + (coords[1] - x.yref)**2)
             for x in potentialMatches]
        d = np.array(d)
        match = potentialMatches[np.argmin(d)]
        return match

    return None


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

        # Get ra/decs of CERTAIN outliers (by distance) for this date
        outliers = GetDistanceOutliers(cluster, date)
        w = GetWCS(cluster, date)

        # Get x- and y-offsets
        if date != baseDate:
            xOffset, yOffset = GetAstrometryOffset(cluster, date, baseDate)
        else:
            xOffset, yOffset = (0, 0)

        for target in data:
            # Compare with other dates, see if target is already found
            xref = binning * target[0] + xOffset
            yref = binning * target[1] + yOffset
            match = MatchTarget(app.cooTol, (xref, yref),
                                BeCandidates + rej_BeCandidates)
            if match is not None:
                # If target refers to an already listed star
                ra = match.ra
                dec = match.dec
                count = match.count
                accept = match in BeCandidates
            else:
                # If target does not match with a listed star
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

            be = BeCandidate()

            be.x = target[0]
            be.y = target[1]
            be.xref = xref
            be.yref = yref
            be.ra = ra
            be.dec = dec
            be.julian = julian
            be.date = date
            be.count = count

            phots = ['Bmag', 'Berr', 'Vmag', 'Verr',
                     'Rmag', 'Rerr', 'Hmag', 'Herr']
            for p, i in zip(phots, range(2, 10)):
                be.phot[p] = target[i]

            be.observed_be = True

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

    # Copy list to stop pass by reference
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
            xOffset, yOffset = (0, 0)

        # Copy list to stop pass by reference
        _targets_to_lookup = list(targets_to_lookup)

        for i in range(1, max([x.count for x in targets_to_lookup]) + 1):
            candidate = [x for x in _targets_to_lookup if x.count == i]
            if any(x.date == date for x in candidate):
                _targets_to_lookup = [x for x in _targets_to_lookup
                                      if x not in candidate]

        for target in data:
            xref = binning * target[0] + xOffset
            yref = binning * target[1] + yOffset
            match = MatchTarget(app.cooTol, (xref, yref),
                                _targets_to_lookup)
            if match is None:
                continue

            be = BeCandidate()

            be.x = target[0]
            be.y = target[1]
            be.xref = xref
            be.yref = yref
            be.ra = match.ra
            be.dec = match.dec
            be.julian = julian
            be.date = date
            be.count = match.count

            phots = ['Bmag', 'Berr', 'Vmag', 'Verr',
                     'Rmag', 'Rerr', 'Hmag', 'Herr']
            for p, i in zip(phots, range(2, 10)):
                be.phot[p] = target[i]

            be.observed_be = False

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
    filename = 'photometry/' + cluster + '/extinctions.dat'
    A_b, A_v, A_r = np.loadtxt(filename).tolist()

    B_VMin = app.B_VMin + A_v - A_b
    B_VMax = app.B_VMax + A_v - A_b

    for target in data:
        b_v = target[2] - target[4]
        if B_VMin < b_v < B_VMax and target[4] <= 13.51:
            BData.append(target)

    # Pick out observed Be stars to result in ONLY B-type stars
    baseDate = ListDates(cluster)[0]
    xOffset, yOffset = GetAstrometryOffset(cluster, date, baseDate)
    binning = Binning(cluster, date)

    for b in reversed(BData):
        xref = b[0] * binning + xOffset
        yref = b[1] * binning + yOffset
        for be in BeCandidates + rej_BeCandidates:
            if abs(xref - be.xref) <= app.cooTol and \
               abs(yref - be.yref) <= app.cooTol:
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
    data = sorted(belist, key=lambda x: (x.count, x.date))

    count = [x.count for x in data]
    for i in range(1, max(count) + 1):
        if i not in count:
            count = [x - 1 for x in count if x > i]
    for i in range(len(data)):
        data[i].count = count[i]

    distance_mean, distance_std = GetRParams(cluster)

    # Open summary files
    if belist == BeCandidates:
        filename = 'output/' + cluster + '/BeList.txt'
    else:
        filename = 'output/' + cluster + '/rejected_BeList.txt'
    list_file = open(filename, 'w')

    for target in data:
        # Get numerical excess
        filename = 'output/' + cluster + '/' + target.date + \
                   '/thresholds.dat'
        thresholds = np.loadtxt(filename)
        if app.threshold_type == 'Constant':
            slope = thresholds[0][0]
            intercept = thresholds[0][1]
        elif app.threshold_type == 'Linear':
            slope = thresholds[1][0]
            intercept = thresholds[1][1]

        r_h = target.phot['Rmag'] - target.phot['Hmag']
        b_v = target.phot['Bmag'] - target.phot['Vmag']
        r_h0 = slope * b_v + intercept
        excess = r_h - r_h0

        # Write to file
        t = '%s-WBBe%d' % (cluster, target.count) + '\t' + \
            '%.3f' % target.GetAbsMag(distance_mean) + '\t' + \
            target.GetSpectralType(distance_mean) + '\t' + \
            '%.10f' % target.ra + '\t' + \
            '%.10f' % target.dec + '\t' + \
            '%.10f' % target.julian + '\t' + \
            '%.3f' % target.phot['Bmag'] + '\t' + \
            '%.3f' % target.phot['Berr'] + '\t' + \
            '%.3f' % target.phot['Vmag'] + '\t' + \
            '%.3f' % target.phot['Verr'] + '\t' + \
            '%.3f' % target.phot['Rmag'] + '\t' + \
            '%.3f' % target.phot['Rerr'] + '\t' + \
            '%.3f' % target.phot['Hmag'] + '\t' + \
            '%.3f' % target.phot['Herr'] + '\t' + \
            '%.3f' % excess + '\n'

        list_file.write(t)

    list_file.close()

    if belist != BeCandidates:
        return

    be_spectralTypes = []
    for i in range(min(count), max(count)):
        m = [x.GetAbsMag(distance_mean) for x in data if x.count == i]
        s = SpectralType(np.mean(m), distance_mean)
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
        b_spectralTypes.append(SpectralType(target[4], distance_mean))

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
    for candidate in (x for x in BeCandidates if x.date == date):
        for target in phot:
            if target[0] == candidate.x and target[1] == candidate.y:
                Be_in.append(target)
                phot.remove(target)
                break

    # Be candidates outside cluster
    Be_out = []
    for candidate in (x for x in rej_BeCandidates if x.date == date):
        for target in phot:
            if target[0] == candidate.x and target[1] == candidate.y:
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

        filename = 'photometry/' + cluster + '/extinctions.dat'
        A_b, A_v, A_r = np.loadtxt(filename).tolist()

        B_VMin = app.B_VMin + A_v - A_b
        B_VMax = app.B_VMax + A_v - A_b

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
            ax.set_ylim([18.5 - A_v, 8.5 - A_v])

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

            # linex = np.array([B_VMin, B_VMax])
            linex = np.array([B_VMin - 0.1, B_VMax + 2.5])
            liney = slope * linex + intercept
            ax.plot(linex, liney, '--', color='#ff5151', label='Be Threshold',
                    linewidth=6)

            output = '2CD_detailed.png'

        # Shared settings
        ax.set_xlabel('B-V')
        ax.set_xlim([B_VMin - 0.1, B_VMax + 2.5])

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
