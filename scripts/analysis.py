import numpy as np
from observations import Observations
from astrometry import Astrometry
from spectral_type import SpectralType
from distance import Distance
import matplotlib.pyplot as plt
# from matplotlib.ticker import MultipleLocator, FormatStrFormatter
from astropy.io import fits
from astropy import wcs
from astropy.coordinates import SkyCoord
from astropy import units as u
import warnings


class Analysis:

    def __init__(self, cluster, app):
        self.cluster = cluster
        self.app = app

        # Open files
        self.File_Summary = open("../output/" + self.cluster + "/summary_" + self.app.phot_type + ".txt", 'w')
        self.BeList = open("../output/" + self.cluster + "/BeList_" + self.app.phot_type + ".txt", 'w')

        # Compile list of Be stars
        self.BeCandidates = []
        self.GetDistanceStats()
        self.CompileBeLists()
        self.ExcludeSingleEntries()
        self.FilterDistances()
        self.FindCorrespondingTargets()

        # Write to files
        self.NightSummary()
        self.BeSummary()

        # Close files
        self.File_Summary.close()
        self.BeList.close()

    def GetDistanceStats(self):
        Distance(self.cluster).process()

        filename = '../output/' + self.cluster + '/distance_params.dat'
        params = np.loadtxt(filename)
        self.distance_mean = params[0]
        self.distance_std = params[1]

    def CompileBeLists(self):
        print("\nCompiling Be list...")

        count = 1
        baseDate = Observations().ListDates(self.cluster)[0]

        for date in Observations().ListDates(self.cluster):
            data = np.loadtxt("../output/" + self.cluster + "/" + date + "/belist_" + self.app.phot_type + ".dat")

            # Get header info
            F = fits.getheader("../photometry/" + self.cluster + "/" + date + "/B1.fits")
            julian = F["JD"]
            binning = F["XBINNING"]
            ra = F['RA']
            dec = F['DEC']
            # Reformat ra/dec
            coo = SkyCoord(ra + dec, unit=(u.hourangle, u.deg))
            ra = coo.ra.deg
            dec = coo.dec.deg

            # Skip process if binning more than 1 because inconsistent (still will be checked in 'FindCorrespondingTargets()' later, however)
            # if binning > 1:
            #     continue

            # Get ra/decs of CERTAIN outliers (by distance) for this date
            distanceCheck = True
            try:
                filename = '../photometry/' + self.cluster + '/' + date + '/phot_dists.csv'
                distanceData = np.genfromtxt(filename, skip_header=1, usecols=(10, 98, 99), delimiter=',')   # parallax, ra, dec

                distanceData = np.array([x for x in distanceData if x[0] > 0])   # don't use negative parallax
                for i in range(0, len(distanceData)):
                    distanceData[i][0] = 1 / distanceData[i][0]   # convert parallax [mas] to distance [kpc]
                distanceData = set(tuple(x) for x in distanceData)
                distanceData = [list(x) for x in distanceData]   # list of unique data

                radec = [[x[1], x[2]] for x in distanceData]
                radec = set(tuple(x) for x in radec)
                radec = [list(x) for x in radec]   # list of unique ra/dec for iteration

                outliers = []

                sigma_coeff = 3
                for target in radec:
                    d = []
                    for line in [x for x in distanceData if x[1] == target[0] and x[2] == target[1]]:
                        print(line)
                        d.append(line[0])
                    if not any(abs(x - self.distance_mean) < sigma_coeff * self.distance_std for x in d):
                        outliers.append(target)
            except IOError:
                print("Note: Data on distances not found for " + date)
                distanceCheck = False

            # Read WCS file for RA and Dec transformations
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                try:
                    w = wcs.WCS("../photometry/" + self.cluster + "/" + date + "/B1_wcs.fits")
                except IOError:
                    print("\nError: Retrieve 'B1_wcs.fits' file for " + self.cluster + " on " + date + " before calculating exact RA and Dec values.  Image center values added as placeholder.")

            if date == baseDate:
                for target in data:
                    try:
                        radec = w.all_pix2world(target[0], target[1], 0)
                        ra = float(radec[0])
                        dec = float(radec[1])
                    except Exception:
                        pass

                    # Exclude if ra/dec corresponds to distance outlier
                    if distanceCheck and [ra, dec] in outliers:
                        continue

                    be = []
                    be.extend([
                        target[0],                       # 0 - x on image
                        target[1],                       # 1 - y on image
                        target[0] * binning,             # 2 - x ref
                        target[1] * binning,             # 3 - y ref
                        "present",                       # 4 - transient status
                        ra,                              # 5 - ra
                        dec,                             # 6 - dec
                        julian,                          # 7 - julian date
                        date,                            # 8 - date
                        count,                           # 9 - count
                        # 10-17 - photometry
                        target[2], target[3], target[4], target[5], target[6], target[7], target[8], target[9],
                        True                             # 18 - be star check
                    ])

                    self.BeCandidates.append(be)
                    count += 1
            else:
                # Get x- and y-offsets
                offsets = Astrometry().GetOffset(self.cluster, date, baseDate)
                xOffset = offsets[0]
                yOffset = offsets[1]

                # Compare with other dates, see if target is already found
                for target in data:
                    be = []
                    xRef = binning * target[0] + xOffset
                    yRef = binning * target[1] + yOffset
                    newCandidate = True
                    for candidate in self.BeCandidates:
                        # if they refer to the same star
                        if abs(candidate[2] - xRef) <= self.app.cooTol and \
                                abs(candidate[3] - yRef) <= self.app.cooTol:
                            # Create Be candidate object
                            be.extend([
                                target[0],                       # 0 - x on image
                                target[1],                       # 1 - y on image
                                xRef,                            # 2 - x ref
                                yRef,                            # 3 - y ref
                                "present",                       # 4 - transient status
                                candidate[5],                    # 5 - ra
                                candidate[6],                    # 6 - dec
                                julian,                          # 7 - julian date
                                date,                            # 8 - date
                                candidate[9],                    # 9 - count
                                # 10-17 - photometry
                                target[2], target[3], target[4], target[5], target[6], target[7], target[8], target[9],
                                True                             # 18 - be star check
                            ])

                            self.BeCandidates.append(be)
                            newCandidate = False
                            break

                    if newCandidate:
                        try:
                            radec = w.all_pix2world(target[0], target[1], 0)
                            ra = float(radec[0])
                            dec = float(radec[1])
                        except Exception:
                            pass

                        # Exclude if ra/dec corresponds to distance outlier
                        if distanceCheck and [ra, dec] in outliers:
                            continue

                        # Create Be candidate object
                        be = []
                        be.extend([
                            target[0],                       # 0 - x on image
                            target[1],                       # 1 - y on image
                            xRef,                            # 2 - x ref
                            yRef,                            # 3 - y ref
                            "present",                       # 4 - transient status
                            ra,                              # 5 - ra
                            dec,                             # 6 - dec
                            julian,                          # 7 - julian date
                            date,                            # 8 - date
                            count,                           # 9 - count
                            # 10-17 - photometry
                            target[2], target[3], target[4], target[5], target[6], target[7], target[8], target[9],
                            True                             # 18 - be star check
                        ])

                        self.BeCandidates.append(be)
                        count += 1

    def ExcludeSingleEntries(self):
        count = [x[9] for x in self.BeCandidates]

        for i in range(0, max(count)):
            if count.count(i) == 1:
                self.BeCandidates = [x for x in self.BeCandidates if x[9] != i]

    def FindCorrespondingTargets(self):
        # Look for non-Be data from previous dates that corresponds to this Be candidate
        baseDate = Observations().ListDates(self.cluster)[0]

        for candidate in self.BeCandidates:
            if candidate[18]:
                for date in Observations().ListDates(self.cluster):
                    be = []
                    # If there aren't any others with same identifier/count and the same date
                    if not any(x[9] == candidate[9] and x[8] == date for x in self.BeCandidates):
                        check_data = np.loadtxt("../output/" + self.cluster + "/" + date + "/phot_" + self.app.phot_type + "_scaled.dat")

                        # Get header info
                        F = fits.getheader("../photometry/" + self.cluster + "/" + date + "/B1.fits")
                        julian = F["JD"]
                        binning = F["XBINNING"]

                        # Get x- and y-offsets
                        if date != baseDate:
                            offsets = Astrometry().GetOffset(self.cluster, date, baseDate)
                            xOffset = offsets[0]
                            yOffset = offsets[1]
                        else:
                            xOffset = 0
                            yOffset = 0

                        for item in check_data:
                            if abs((binning * item[0] + xOffset) - candidate[2]) <= self.app.cooTol and \
                                    abs((binning * item[1] + yOffset) - candidate[3]) <= self.app.cooTol:
                                be.extend([
                                    item[0],                         # 0 - x on image
                                    item[1],                         # 1 - y on image
                                    binning * item[0] + xOffset,     # 2 - x ref
                                    binning * item[1] + yOffset,     # 3 - y ref
                                    "absent",                        # 4 - transient status
                                    candidate[5],                    # 5 - ra
                                    candidate[6],                    # 6 - dec
                                    julian,                          # 7 - julian date
                                    date,                            # 8 - date
                                    candidate[9],                    # 9 - count
                                    # 10-17 - photometry
                                    item[2], item[3], item[4], item[5], item[6], item[7], item[8], item[9],
                                    False                            # 18 - be star check
                                ])
                                self.BeCandidates.append(be)
                                break

    def FindBStars(self, date):
        print("Finding B-type stars for " + date)

        path = "../output/" + self.cluster + "/" + date + "/"
        data = np.loadtxt(path + "phot_" + self.app.phot_type + "_scaled.dat")
        BData = []

        # Find all B- and Be-type stars
        for target in data:
            b_v = target[2] - target[4]
            if b_v < self.app.B_VMax and b_v > self.app.B_VMin and target[4] <= 13.51:
                BData.append(target)

        # Pick out observed Be stars to result in ONLY B-type stars
        baseDate = Observations().ListDates(self.cluster)[0]
        offsets = Astrometry().GetOffset(self.cluster, date, baseDate)
        xOffset = offsets[0]
        yOffset = offsets[1]
        F = fits.getheader("../photometry/" + self.cluster + "/" + date + "/B1.fits")
        binning = F["XBINNING"]

        for b in BData:
            xref = b[0] * binning + xOffset
            yref = b[1] * binning + yOffset
            for be in self.BeCandidates:
                if abs(xref - be[2]) <= self.app.cooTol and abs(yref - be[3]) <= self.app.cooTol:
                    np.delete(BData, np.argwhere(BData == b))

        return BData

    def BeValues(self, date):
        # Read data
        beData = np.loadtxt("../output/" + self.cluster + "/" + date + "/beList_" + self.app.phot_type + ".dat")
        NumBe = len(beData)

        bData = self.FindBStars(date)
        NumB = len(bData)

        Be_ratio = NumBe / (NumB + NumBe)

        self.File_Summary.write("Be candidates:                " + str(NumBe) + "\n")
        self.File_Summary.write("B stars:                      " + str(NumB) + "\n")
        self.File_Summary.write("Be ratio:                     " + "%.3f" % Be_ratio + "\n")
        self.File_Summary.write("\n\n")

        if NumB >= self.mostB:
            self.mostB = NumB
            self.mostB_date = date

    def NightSummary(self):
        self.File_Summary.write("================================================\n")
        self.File_Summary.write("                " + self.cluster + " Summary                  \n")
        self.File_Summary.write("================================================\n\n\n")

        self.mostB = 0
        dates = Observations().ListDates(self.cluster)
        for date in dates:
            self.File_Summary.write("                   " + date + "                   \n")
            self.File_Summary.write("------------------------------------------------\n")

            # Technical information
            filename = "../photometry/" + self.cluster + "/" + date + "/B1.fits"
            G = fits.getheader(filename)
            self.File_Summary.write("Binning: " + str(G["XBINNING"]) + '\n')

            self.File_Summary.write("Exposures: ")
            files = ['B1', 'B3', 'V1', 'V3', 'R1', 'R3', 'H1', 'H3']
            for x in files:
                G = fits.getheader("../photometry/" + self.cluster + "/" + date + "/" + x + ".fits")
                self.File_Summary.write(x + ": " + '%d' % G["EXPTIME"])
                if x != files[len(files) - 1]:
                    self.File_Summary.write(" | ")

            self.File_Summary.write('\n\n')

            # Threshold information
            filename = "../output/" + self.cluster + "/" + date + "/thresholds_" + self.app.phot_type + ".dat"
            thresholds = np.loadtxt(filename)
            self.File_Summary.write("Thresholds:\n")
            self.File_Summary.write("   Constant: " + "%.3f" % thresholds[0][1] + "     Linear: " + "%.3f" % thresholds[1][0] + ", " + "%.3f" % thresholds[1][1] + "\n\n")

            # Be candidate information
            self.BeValues(date)

    def BeSummary(self):
        data = sorted(self.BeCandidates, key=lambda x: (x[9], x[8]))   # Sort by identifier (count) then date

        ra = [x[5] for x in data]
        dec = [x[6] for x in data]
        julian = [x[7] for x in data]
        date = [x[8] for x in data]

        count = [x[9] for x in data]
        for i in range(0, max(count)):
            if i not in count:
                for j in range(0, len(count)):
                    if count[j] > i:
                        count[j] -= 1

        bmag = [x[10] for x in data]
        berr = [x[11] for x in data]
        vmag = [x[12] for x in data]
        verr = [x[13] for x in data]
        rmag = [x[14] for x in data]
        rerr = [x[15] for x in data]
        hmag = [x[16] for x in data]
        herr = [x[17] for x in data]

        absVmag = [SpectralType(self.distance_mean).AbsMag(x) for x in vmag]   # 19
        spectralTypes = [SpectralType(self.distance_mean).GetSpectralType(x) for x in vmag]   # 20
        for i in range(0, len(data)):
            data[i].extend([absVmag[i], spectralTypes[i]])

        for i in range(0, len(data)):
            # Get numerical excess
            try:
                filename = "../output/" + self.cluster + "/" + date[i] + "/thresholds_" + self.app.phot_type + ".dat"
                thresholds = np.loadtxt(filename)
                if self.app.threshold_type == "Constant":
                    slope = thresholds[0][0]
                    intercept = thresholds[0][1]
                elif self.app.threshold_type == "Linear":
                    slope = thresholds[1][0]
                    intercept = thresholds[1][1]
            except Exception:
                print("Error: Threshold calculations are required.")

            r_h = rmag[i] - hmag[i]
            b_v = bmag[i] - vmag[i]
            r_h0 = slope * b_v + intercept
            excess = r_h - r_h0

            # Write to file
            self.BeList.write(self.cluster + "-WBBe" + str(count[i]) + "\t" +
                              "%.3f" % absVmag[i] + "\t" +
                              spectralTypes[i] + "\t" +
                              "%.10f" % ra[i] + "\t" +
                              "%.10f" % dec[i] + "\t" +
                              "%.10f" % julian[i] + "\t" +
                              "%.3f" % bmag[i] + "\t" +
                              "%.3f" % berr[i] + "\t" +
                              "%.3f" % vmag[i] + "\t" +
                              "%.3f" % verr[i] + "\t" +
                              "%.3f" % rmag[i] + "\t" +
                              "%.3f" % rerr[i] + "\t" +
                              "%.3f" % hmag[i] + "\t" +
                              "%.3f" % herr[i] + "\t" +
                              "%.3f" % excess + "\n")

        be_spectralTypes = []
        for i in range(0, max(count)):
            m = [x[12] for x in data if x[9] == i]
            be_spectralTypes.append(SpectralType(self.distance_mean).GetSpectralType(np.mean(m)))

        be_type_unknown = [x for x in be_spectralTypes if x == "--"]
        be_type_O = [x for x in be_spectralTypes if x[0] == "O"]
        be_type_B0_B3 = [x for x in be_spectralTypes if x == "B0" or x == "B1" or x == "B2" or x == "B3"]
        be_type_B4_B5 = [x for x in be_spectralTypes if x == "B4" or x == "B5"]
        be_type_B6_B9 = [x for x in be_spectralTypes if x == "B6" or x == "B7" or x == "B8" or x == "B9"]
        be_type_A = [x for x in be_spectralTypes if x[0] == "A"]

        frequencies = [len(be_type_unknown), len(be_type_O), len(be_type_B0_B3), len(be_type_B4_B5), len(be_type_B6_B9), len(be_type_A)]
        names = ["<O6", "O6-O8", "B0-B3", "B4-B5", "B6-B9", "A"]

        # Plot spectral type histogram

        plt.figure(figsize=(12, 9))

        x_coordinates = np.arange(len(frequencies))
        plt.bar(x_coordinates, frequencies, align='center', color="#3f3f3f")
        # plt.title(filter + " Magnitude Scaling Differences")
        plt.xlabel("Spectral Type", fontsize=36)
        plt.ylabel("Frequency", fontsize=36)

        plt.axes().xaxis.set_major_locator(plt.FixedLocator(x_coordinates))
        plt.axes().xaxis.set_major_formatter(plt.FixedFormatter(names))

        plt.axes().tick_params('both', length=12, width=4, which='major', top=True, right=True, direction='in', pad=6, labelsize=30)
        plt.axes().tick_params('both', length=8, width=3, which='minor', top=True, right=True, direction='in')

        plt.axes().spines['top'].set_linewidth(4)
        plt.axes().spines['right'].set_linewidth(4)
        plt.axes().spines['bottom'].set_linewidth(4)
        plt.axes().spines['left'].set_linewidth(4)

        plt.tight_layout()

        filename = "../output/" + self.cluster + "/spectral_types_" + self.app.phot_type + ".png"
        plt.savefig(filename)
        plt.clf()

        # Be ratios by spectral type
        BData = self.FindBStars(self.mostB_date)
        b_spectralTypes = []
        for i in range(0, len(BData)):
            b_spectralTypes.append(SpectralType(self.distance_mean).GetSpectralType(BData[i][4]))

        b_type_unknown = len([x for x in b_spectralTypes if x[0] == "--"])
        b_type_O = len([x for x in b_spectralTypes if x[0] == "O"])
        b_type_B0_B3 = len([x for x in b_spectralTypes if x == "B0" or x == "B1" or x == "B2" or x == "B3"])
        b_type_B4_B5 = len([x for x in b_spectralTypes if x == "B4" or x == "B5"])
        b_type_B6_B9 = len([x for x in b_spectralTypes if x == "B6" or x == "B7" or x == "B8" or x == "B9"])
        b_type_A = len([x for x in b_spectralTypes if x[0] == "A"])

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

        filename = "../output/" + self.cluster + "/ratios_by_spec_type.txt"
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
