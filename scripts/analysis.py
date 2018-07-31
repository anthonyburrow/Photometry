import numpy as np
from observations import Observations
from astrometry import Astrometry
from astropy.io import fits
from astropy import wcs
import warnings


class Analysis:

    def __init__(self, cluster, app):
        self.cluster = cluster
        self.app = app
        self.BeCandidates = []

    def Values(self, date):
        # Read data
        data = np.loadtxt("../output/" + self.cluster + "/" + date + "/phot_" + self.app.phot_type + ".dat")
        beData = np.loadtxt("../output/" + self.cluster + "/" + date + "/beList_" + self.app.phot_type + ".dat")

        data_lowError = np.loadtxt("../output/" + self.cluster + "/" + date + "/phot_" + self.app.phot_type + "_lowError.dat")
        beData_lowError = np.loadtxt("../output/" + self.cluster + "/" + date + "/beList_" + self.app.phot_type + "_lowError.dat")

        # Number of Be-type stars
        self.NumBe = len(beData)
        self.NumBe_lowError = len(beData_lowError)

        # Number of B-type stars
        BTotal = []
        BTotal_lowError = []

        for target in data:
            B_V = target[2] - target[4]
            if B_V > self.app.B_VMin and B_V < self.app.B_VMax and \
                    target[4] <= 13.51:
                BTotal.append(target)
        for target in data_lowError:
            B_V = target[2] - target[4]
            if B_V > self.app.B_VMin and B_V < self.app.B_VMax and \
                    target[4] <= 13.51:
                BTotal_lowError.append(target)

        self.NumBTotal = len(BTotal)
        self.NumBTotal_lowError = len(BTotal_lowError)

        # Ratio of Be to Be + B (total)
        self.Be_ratio = self.NumBe / self.NumBTotal
        self.Be_ratio_lowError = self.NumBe_lowError / self.NumBTotal_lowError

    def CompileBeLists(self):
        count = 1

        baseDate = Observations().ListDates(self.cluster)[0]
        for date in Observations().ListDates(self.cluster):
            data = np.loadtxt("../output/" + self.cluster + "/" + date + "/belist_" + self.app.phot_type + ".dat")

            # Get header info
            F = fits.getheader("../photometry/" + self.cluster + "/" + date + "/B1.fits")
            julian = F["JD"]
            binning = F["XBINNING"]
            ra = F["CRVAL1"]
            dec = F["CRVAL1"]

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
                        ra, dec = w.all_pix2world(target[0], target[1], 0)
                    except Exception:
                        pass

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
                # Get x- and y- offsets
                filename = "../output/" + self.cluster + "/" + date + "/astrometry_offsets.txt"
                try:
                    offsets = np.loadtxt(filename)
                except IOError:
                    offsets = Astrometry().GetOffset(self.cluster, date, baseDate)
                    with open(filename, 'w') as F:
                        np.savetxt(F, offsets, fmt="%.3f")

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
                            # Determine if it was gained on this date
                            # prev_date_ind = Observations().ListDates(self.cluster).index(date) - 1
                            # prev_date = Observations().ListDates(self.cluster)[prev_date_ind]
                            # for item in self.BeCandidates:
                            #     if item[9] == candidate[9] and item[8] == prev_date and item[18]:
                            #         transient = "  --  "
                            #         break

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
                            ra, dec = w.all_pix2world(target[0], target[1], 0)
                        except Exception:
                            pass

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

    def FindCorrespondingTargets(self):
        # Look for non-Be data from previous dates that corresponds to this Be candidate
        baseDate = Observations().ListDates(self.cluster)[0]

        for candidate in self.BeCandidates:
            if candidate[18]:
                for date in Observations().ListDates(self.cluster):
                    be = []
                    # If there aren't any others with same identifier/count and the same date
                    if not any(x[9] == candidate[9] and x[8] == date for x in self.BeCandidates):
                        check_data = np.loadtxt("../output/" + self.cluster + "/" + date + "/phot_" + self.app.phot_type + ".dat")

                        # Get header info
                        F = fits.getheader("../photometry/" + self.cluster + "/" + date + "/B1.fits")
                        julian = F["JD"]
                        binning = F["XBINNING"]

                        # Get x- and y- offsets
                        if date != baseDate:
                            filename = "../output/" + self.cluster + "/" + date + "/astrometry_offsets.txt"
                            offsets = np.loadtxt(filename)
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

    def Summary(self):
        with open("../output/" + self.cluster + "/summary_" + self.app.phot_type + ".txt", 'w') as F:
            F.write("================================================\n")
            F.write("                " + self.cluster + " Summary                  \n")
            F.write("================================================\n\n\n")

            for date in Observations().ListDates(self.cluster):
                F.write("                   " + date + "                   \n")
                F.write("------------------------------------------------\n")

                # Technical information
                filename = "../photometry/" + self.cluster + "/" + date + "/B1.fits"
                G = fits.getheader(filename)
                F.write("Binning: " + str(G["XBINNING"]) + '\n')

                F.write("Exposures: ")
                files = ['B1', 'B3', 'V1', 'V3', 'R1', 'R3', 'H1', 'H3']
                for x in files:
                    G = fits.getheader("../photometry/" + self.cluster + "/" + date + "/" + x + ".fits")
                    F.write(x + ": " + '%d' % G["EXPTIME"])
                    if x != files[len(files) - 1]:
                        F.write(" | ")

                F.write('\n\n')

                # Threshold information
                filename = "../output/" + self.cluster + "/" + date + "/thresholds_" + self.app.phot_type + ".dat"
                thresholds = np.loadtxt(filename)
                F.write("Thresholds:\n")
                F.write("   Constant: " + "%.3f" % thresholds[0][1] + "     Linear: " + "%.3f" % thresholds[1][0] + ", " + "%.3f" % thresholds[1][1] + "\n\n")

                # Be candidate information
                self.Values(date)

                F.write("Be candidates:                " + str(self.NumBe) + "\n")
                F.write("Be candidates (low error):    " + str(self.NumBe_lowError) + "\n")
                F.write("B stars:                      " + str(self.NumBTotal - self.NumBe) + "\n")
                F.write("B stars (low error):          " + str(self.NumBTotal_lowError - self.NumBe_lowError) + "\n")
                F.write("Be ratio:                     " + "%.3f" % self.Be_ratio + "\n")
                F.write("Be ratio (low error):         " + "%.3f" % self.Be_ratio_lowError + "\n")
                F.write("\n\n")

        # Write compiled Be list
        self.CompileBeLists()
        self.FindCorrespondingTargets()

        data = self.BeCandidates
        data = sorted(data, key=lambda x: (x[9], x[8]))   # Sort by identifier (count) then date

        transient = [x[4] for x in data]
        ra = [x[5] for x in data]
        dec = [x[6] for x in data]
        julian = [x[7] for x in data]
        date = [x[8] for x in data]
        count = [x[9] for x in data]

        bmag = [x[10] for x in data]
        berr = [x[11] for x in data]
        vmag = [x[12] for x in data]
        verr = [x[13] for x in data]
        rmag = [x[14] for x in data]
        rerr = [x[15] for x in data]
        hmag = [x[16] for x in data]
        herr = [x[17] for x in data]

        filename = "../output/" + self.cluster + "/BeList_" + self.app.phot_type + ".txt"
        with open(filename, 'w') as F:
            for i in range(0, len(data)):
                F.write(self.cluster + "-WBBe" + str(count[i]) + "\t" +
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
                        transient[i] + "\n")

        filename = "../output/" + self.cluster + "/HaExcesses_" + self.app.phot_type + ".txt"
        with open(filename, 'w') as F:
            for i in range(0, len(data)):
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

                F.write(self.cluster + "-WBBe" + str(count[i]) + "\t" +
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
