from .observations import ListDates
from .astrometry import GetAstrometryOffset
from .spectral_type import AbsMag, GetSpectralType
from .distance import ProcessDistances

import matplotlib.pyplot as plt
# from matplotlib.ticker import MultipleLocator, FormatStrFormatter

import numpy as np
from astropy.io import fits
from astropy import wcs
from astropy.coordinates import SkyCoord
from astropy import units as u
import warnings


def ProcessAnalysis(cluster, app):
    ProcessDistances(cluster)

    # Compile list of Be stars
    BeCandidates = []
    rej_BeCandidates = []
    CompileBeLists(cluster, app, BeCandidates, rej_BeCandidates)
    FindCorrespondingTargets(cluster, app, BeCandidates)
    FindCorrespondingTargets(cluster, app, rej_BeCandidates)

    # Write to files
    mostB_date = NightSummary(cluster, app, BeCandidates)
    BeSummary(cluster, app, BeCandidates, mostB_date, BeCandidates)
    BeSummary(cluster, app, rej_BeCandidates, mostB_date, BeCandidates)


def GetDistanceParams(cluster):
    filename = 'output/' + cluster + '/distance_params.dat'
    params = np.loadtxt(filename)

    return params.tolist()


def GetDistanceOutliers(cluster, date):
    outliers = []
    try:
        filename = 'photometry/' + cluster + '/' + date + '/phot_dists.csv'
        distanceData = np.genfromtxt(filename, skip_header=1, usecols=(10, 98, 99), delimiter=',')   # parallax, ra, dec

        distanceData = np.array([x for x in distanceData if x[0] > 0])   # don't use negative parallax
        for target in distanceData:
            target[0] = 1 / target[0]   # convert parallax [mas] to distance [kpc]
        distanceData = set(tuple(x) for x in distanceData)
        distanceData = [list(x) for x in distanceData]   # list of unique data

        radec = [[x[1], x[2]] for x in distanceData]
        radec = set(tuple(x) for x in radec)
        radec = [list(x) for x in radec]   # list of unique ra/dec for iteration

        sigma_coeff = 3
        for target in radec:
            d = []
            for line in [x for x in distanceData if x[1] == target[0] and x[2] == target[1]]:
                d.append(line[0])
            distance_mean, distance_std = GetDistanceParams(cluster)
            if not any(abs(x - distance_mean) < sigma_coeff * distance_std for x in d):
                outliers.append(target)
    except IOError:
        print("Note: Data on distances not found for " + date)

    return outliers


def CompileBeLists(cluster, app, BeCandidates, rej_BeCandidates):
    print("\nCompiling Be list...")

    acc_count = 1
    rej_count = 1
    baseDate = ListDates(cluster)[0]

    for date in ListDates(cluster):
        filename = 'output/' + cluster + '/' + date + '/belist_scaled_' + app.phot_type + '.dat'
        data = np.loadtxt(filename)

        # Get header info
        filename = 'photometry/' + cluster + '/' + date + '/B1.fits'
        F = fits.getheader(filename)
        julian = F['JD']
        binning = F['XBINNING']
        ra = F['RA']
        dec = F['DEC']
        # Reformat ra/dec
        coo = SkyCoord(ra + dec, unit=(u.hourangle, u.deg))
        ra = coo.ra.deg
        dec = coo.dec.deg

        # Get ra/decs of CERTAIN outliers (by distance) for this date
        outliers = GetDistanceOutliers(cluster, date)

        # Read WCS file for RA and Dec transformations
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            try:
                filename = 'photometry/' + cluster + '/' + date + '/B1_wcs.fits'
                w = wcs.WCS(filename)
            except IOError:
                print("\nError: Retrieve 'B1_wcs.fits' file for " + cluster + " on " + date +
                      " before calculating exact RA and Dec values. \
                        Image center values added as placeholder.")

        if date == baseDate:
            for target in data:
                try:
                    radec = w.all_pix2world(target[0], target[1], 0)
                    ra = float(radec[0])
                    dec = float(radec[1])
                except Exception:
                    pass

                if [ra, dec] not in outliers:
                    count = acc_count
                else:
                    count = rej_count

                be = []
                be.extend([
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
                ])

                if [ra, dec] not in outliers:
                    BeCandidates.append(be)
                    acc_count += 1
                else:
                    rej_BeCandidates.append(be)
                    rej_count += 1
        else:
            # Get x- and y-offsets
            xOffset, yOffset = GetAstrometryOffset(cluster, date, baseDate)

            # Compare with other dates, see if target is already found
            for target in data:
                be = []
                xRef = binning * target[0] + xOffset
                yRef = binning * target[1] + yOffset
                for candidate in BeCandidates + rej_BeCandidates:
                    # if they refer to the same star
                    if abs(candidate[2] - xRef) <= app.cooTol and \
                            abs(candidate[3] - yRef) <= app.cooTol:
                        be.extend([
                            target[0],                       # 0 - x on image
                            target[1],                       # 1 - y on image
                            xRef,                            # 2 - x ref
                            yRef,                            # 3 - y ref
                            'present',                       # 4 - transient status
                            candidate[5],                    # 5 - ra
                            candidate[6],                    # 6 - dec
                            julian,                          # 7 - julian date
                            date,                            # 8 - date
                            candidate[9],                    # 9 - count
                            # 10-17 - photometry
                            target[2], target[3], target[4], target[5],
                            target[6], target[7], target[8], target[9],
                            True                             # 18 - be star check
                        ])

                        if candidate in BeCandidates:
                            BeCandidates.append(be)
                        else:
                            rej_BeCandidates.append(be)

                        break
                else:
                    try:
                        radec = w.all_pix2world(target[0], target[1], 0)
                        ra = float(radec[0])
                        dec = float(radec[1])
                    except Exception:
                        pass

                    if [ra, dec] not in outliers:
                        count = acc_count
                    else:
                        count = rej_count

                    # Create Be candidate object
                    be = []
                    be.extend([
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
                    ])

                    if [ra, dec] not in outliers:
                        BeCandidates.append(be)
                        acc_count += 1
                    else:
                        rej_BeCandidates.append(be)
                        rej_count += 1


def FindCorrespondingTargets(cluster, app, belist):
    # Look for non-Be data from previous dates that corresponds to this Be candidate
    baseDate = ListDates(cluster)[0]

    for candidate in [x for x in belist if x[18]]:
        for date in ListDates(cluster):
            be = []
            # If there aren't any others with same identifier/count and the same date
            if not any(x[9] == candidate[9] and x[8] == date for x in belist):
                filename = 'output/' + cluster + '/' + date + '/phot_scaled_' + app.phot_type + '.dat'
                check_data = np.loadtxt(filename)

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

                for item in check_data:
                    if abs((binning * item[0] + xOffset) - candidate[2]) <= app.cooTol and \
                            abs((binning * item[1] + yOffset) - candidate[3]) <= app.cooTol:
                        be.extend([
                            item[0],                         # 0 - x on image
                            item[1],                         # 1 - y on image
                            binning * item[0] + xOffset,     # 2 - x ref
                            binning * item[1] + yOffset,     # 3 - y ref
                            'absent',                        # 4 - transient status
                            candidate[5],                    # 5 - ra
                            candidate[6],                    # 6 - dec
                            julian,                          # 7 - julian date
                            date,                            # 8 - date
                            candidate[9],                    # 9 - count
                            # 10-17 - photometry
                            item[2], item[3], item[4], item[5], item[6], item[7], item[8], item[9],
                            False                            # 18 - be star check
                        ])
                        belist.append(be)
                        break


def FindBStars(cluster, app, date, BeCandidates):
    print("Finding B-type stars for " + date)

    filename = 'output/' + cluster + '/' + date + '/phot_scaled_' + app.phot_type + '.dat'
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
        for be in BeCandidates:
            if abs(xref - be[2]) <= app.cooTol and abs(yref - be[3]) <= app.cooTol:
                BData.remove(b)
                break

    # Get rid of stars not in cluster (by distance)
    outliers = GetDistanceOutliers(cluster, date)

    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        try:
            w = wcs.WCS('photometry/' + cluster + '/' + date + '/B1_wcs.fits')
        except IOError:
            print("\nError: Retrieve 'B1_wcs.fits' file for " + cluster + " on " + date +
                  " before calculating exact RA and Dec values. \
                    Image center values added as placeholder.")

    for target in BData:
        radec = w.all_pix2world(target[0], target[1], 0)
        ra = float(radec[0])
        dec = float(radec[1])

        # If ra/dec corresponds to distance outlier, skip target
        if [ra, dec] in reversed(outliers):
            BData.remove(target)

    return BData


def BeValues(cluster, app, date, summary_file, mostB, BeCandidates):
    # Read data
    filename = 'output/' + cluster + '/' + date + '/beList_scaled_' + app.phot_type + '.dat'
    beData = np.loadtxt(filename)
    NumBe = len(beData)

    bData = FindBStars(cluster, app, date, BeCandidates)
    NumB = len(bData)

    Be_ratio = NumBe / (NumB + NumBe)

    summary_file.write("Be candidates:                " + str(NumBe) + "\n")
    summary_file.write("B stars:                      " + str(NumB) + "\n")
    summary_file.write("Be ratio:                     " + "%.3f" % Be_ratio + "\n")
    summary_file.write("\n\n")

    mostB_new = mostB
    if NumB > mostB[0]:
        mostB_new = [NumB, date]

    return mostB_new


def NightSummary(cluster, app, BeCandidates):
    filename = 'output/' + cluster + '/summary_' + app.phot_type + '.txt'
    summary_file = open(filename, 'w')

    summary_file.write("================================================\n")
    summary_file.write("                " + cluster + " Summary                  \n")
    summary_file.write("================================================\n\n\n")

    mostB = [0, '']
    dates = ListDates(cluster)
    for date in dates:
        summary_file.write("                   " + date + "                   \n")
        summary_file.write("------------------------------------------------\n")

        # Technical information
        filename = 'photometry/' + cluster + '/' + date + '/B1.fits'
        G = fits.getheader(filename)
        summary_file.write("Binning: " + str(G['XBINNING']) + '\n')

        summary_file.write("Exposures: ")
        files = ['B1', 'B3', 'V1', 'V3', 'R1', 'R3', 'H1', 'H3']
        for x in files:
            filename = 'photometry/' + cluster + '/' + date + '/' + x + '.fits'
            G = fits.getheader(filename)
            summary_file.write(x + ": " + '%d' % G['EXPTIME'])
            if x != files[len(files) - 1]:
                summary_file.write(" | ")

        summary_file.write('\n\n')

        # Threshold information
        filename = 'output/' + cluster + '/' + date + '/thresholds_' + app.phot_type + '.dat'
        thresholds = np.loadtxt(filename)
        summary_file.write("Thresholds:\n")
        summary_file.write("   Constant: " + "%.3f" % thresholds[0][1] +
                           "     Linear: " + "%.3f" % thresholds[1][0] + ", " +
                                             "%.3f" % thresholds[1][1] + "\n\n")

        # Be candidate information
        mostB = BeValues(cluster, app, date, summary_file, mostB, BeCandidates)

    summary_file.close()

    return mostB[1]


def BeSummary(cluster, app, belist, mostB_date, BeCandidates):
    if not belist:
        return

    data = sorted(belist, key=lambda x: (x[9], x[8]))   # Sort by identifier (count) then date

    ra = [x[5] for x in data]
    dec = [x[6] for x in data]
    julian = [x[7] for x in data]
    date = [x[8] for x in data]

    count = [x[9] for x in data]
    for i in range(1, max(count)):
        if i not in count:
            for j in count:
                if j > i:
                    j -= 1

    bmag = [x[10] for x in data]
    berr = [x[11] for x in data]
    vmag = [x[12] for x in data]
    verr = [x[13] for x in data]
    rmag = [x[14] for x in data]
    rerr = [x[15] for x in data]
    hmag = [x[16] for x in data]
    herr = [x[17] for x in data]

    distance_mean, distance_std = GetDistanceParams(cluster)

    absVmag = [AbsMag(x, distance_mean) for x in vmag]   # 19
    spectralTypes = [GetSpectralType(x, distance_mean) for x in vmag]   # 20
    for i in range(0, len(data)):
        data[i].extend([absVmag[i], spectralTypes[i]])

    # Open summary files
    if belist == BeCandidates:
        filename = 'output/' + cluster + '/BeList_' + app.phot_type + '.txt'
    else:
        filename = 'output/' + cluster + '/rejected_BeList_' + app.phot_type + '.txt'
    list_file = open(filename, 'w')

    for i in range(0, len(data)):
        # Get numerical excess
        try:
            filename = 'output/' + cluster + '/' + date[i] + '/thresholds_' + app.phot_type + '.dat'
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
        wrt = cluster + '-WBBe' + str(count[i]) + '\t' + \
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

        list_file.write(wrt)

    list_file.close()

    if belist != BeCandidates:
        return

    be_spectralTypes = []
    for i in range(0, max(count)):
        m = [x[12] for x in data if x[9] == i]
        s = GetSpectralType(np.mean(m), distance_mean)
        be_spectralTypes.append(s)

    be_type_unknown = [x for x in be_spectralTypes if x == '--']
    be_type_O = [x for x in be_spectralTypes if x[0] == 'O']
    be_type_B0_B3 = [x for x in be_spectralTypes if x == 'B0' or x == 'B1' or x == 'B2' or x == 'B3']
    be_type_B4_B5 = [x for x in be_spectralTypes if x == 'B4' or x == 'B5']
    be_type_B6_B9 = [x for x in be_spectralTypes if x == 'B6' or x == 'B7' or x == 'B8' or x == 'B9']
    be_type_A = [x for x in be_spectralTypes if x[0] == 'A']

    frequencies = [len(be_type_unknown), len(be_type_O), len(be_type_B0_B3), len(be_type_B4_B5), len(be_type_B6_B9), len(be_type_A)]
    names = ['<O6', 'O6-O8', 'B0-B3', 'B4-B5', 'B6-B9', 'A']

    # Plot spectral type histogram

    plt.figure(figsize=(12, 9))

    x_coordinates = np.arange(len(frequencies))
    plt.bar(x_coordinates, frequencies, align='center', color='#3f3f3f')
    # plt.title(filter + " Magnitude Scaling Differences")
    plt.xlabel('Spectral Type', fontsize=36)
    plt.ylabel('Frequency', fontsize=36)

    plt.axes().xaxis.set_major_locator(plt.FixedLocator(x_coordinates))
    plt.axes().xaxis.set_major_formatter(plt.FixedFormatter(names))

    plt.axes().tick_params('both', length=12, width=4, which='major', top=True, right=True, direction='in', pad=6, labelsize=30)
    plt.axes().tick_params('both', length=8, width=3, which='minor', top=True, right=True, direction='in')

    plt.axes().spines['top'].set_linewidth(4)
    plt.axes().spines['right'].set_linewidth(4)
    plt.axes().spines['bottom'].set_linewidth(4)
    plt.axes().spines['left'].set_linewidth(4)

    plt.tight_layout()

    filename = 'output/' + cluster + '/spectral_types_' + app.phot_type + '.png'
    plt.savefig(filename)
    plt.clf()

    # Be ratios by spectral type
    BData = FindBStars(cluster, app, mostB_date, BeCandidates)
    b_spectralTypes = []
    for target in BData:
        b_spectralTypes.append(GetSpectralType(target[4], distance_mean))

    b_type_unknown = len([x for x in b_spectralTypes if x[0] == '--'])
    b_type_O = len([x for x in b_spectralTypes if x[0] == 'O'])
    b_type_B0_B3 = len([x for x in b_spectralTypes if x == 'B0' or x == 'B1' or x == 'B2' or x == 'B3'])
    b_type_B4_B5 = len([x for x in b_spectralTypes if x == 'B4' or x == 'B5'])
    b_type_B6_B9 = len([x for x in b_spectralTypes if x == 'B6' or x == 'B7' or x == 'B8' or x == 'B9'])
    b_type_A = len([x for x in b_spectralTypes if x[0] == 'A'])

    try:
        ratio_unknown = len(be_type_unknown) / (len(be_type_unknown) + b_type_unknown)
    except ZeroDivisionError:
        ratio_unknown = 0
    try:
        ratio_O = len(be_type_O) / (len(be_type_O) + b_type_O)
    except ZeroDivisionError:
        ratio_O = 0
    try:
        ratio_B0_B3 = len(be_type_B0_B3) / (len(be_type_B0_B3) + b_type_B0_B3)
    except ZeroDivisionError:
        ratio_B0_B3 = 0
    try:
        ratio_B4_B5 = len(be_type_B4_B5) / (len(be_type_B4_B5) + b_type_B4_B5)
    except ZeroDivisionError:
        ratio_B4_B5 = 0
    try:
        ratio_B6_B9 = len(be_type_B6_B9) / (len(be_type_B6_B9) + b_type_B6_B9)
    except ZeroDivisionError:
        ratio_B6_B9 = 0
    try:
        ratio_A = len(be_type_A) / (len(be_type_A) + b_type_A)
    except ZeroDivisionError:
        ratio_A = 0

    filename = 'output/' + cluster + '/ratios_by_spec_type.txt'
    with open(filename, 'w') as F:
        F.write("Unknown:\n")
        F.write("    Be: " + str(len(be_type_unknown)) + '\n')
        F.write("    B: " + str(b_type_unknown) + '\n')
        F.write("    Ratio: " + '%.3f' % ratio_unknown + '\n\n')
        F.write("Type O:\n")
        F.write("    Be: " + str(len(be_type_O)) + '\n')
        F.write("    B: " + str(b_type_O) + '\n')
        F.write("    Ratio: " + '%.3f' % ratio_O + '\n\n')
        F.write("Type B0-B3:\n")
        F.write("    Be: " + str(len(be_type_B0_B3)) + '\n')
        F.write("    B: " + str(b_type_B0_B3) + '\n')
        F.write("    Ratio: " + '%.3f' % ratio_B0_B3 + '\n\n')
        F.write("Type B4-B5:\n")
        F.write("    Be: " + str(len(be_type_B4_B5)) + '\n')
        F.write("    B: " + str(b_type_B4_B5) + '\n')
        F.write("    Ratio: " + '%.3f' % ratio_B4_B5 + '\n\n')
        F.write("Type B6-B9:\n")
        F.write("    Be: " + str(len(be_type_B6_B9)) + '\n')
        F.write("    B: " + str(b_type_B6_B9) + '\n')
        F.write("    Ratio: " + '%.3f' % ratio_B6_B9 + '\n\n')
        F.write("Type A:\n")
        F.write("    Be: " + str(len(be_type_A)) + '\n')
        F.write("    B: " + str(b_type_A) + '\n')
        F.write("    Ratio: " + '%.3f' % ratio_A + '\n\n')
