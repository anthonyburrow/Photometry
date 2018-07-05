import numpy as np
from observations import Observations
from astrometry import Astrometry
from astropy.io import fits


class Analysis:

    def __init__(self, cluster, app):
        self.cluster = cluster
        self.app = app

    def Values(self, date):
        # Read data
        data = np.loadtxt("../output/" + self.cluster + "/" + date + "/phot_" + self.app.phot_type + ".dat")
        beData = np.loadtxt("../output/" + self.cluster + "/" + date + "/beList.dat")

        data_lowError = np.loadtxt("../output/" + self.cluster + "/" + date + "/phot_" + self.app.phot_type + "_lowError.dat")
        beData_lowError = np.loadtxt("../output/" + self.cluster + "/" + date + "/beList_lowError.dat")

        # Number of Be-type stars
        self.NumBe = len(beData)
        self.NumBe_lowError = len(beData_lowError)

        # Number of B-type stars
        BTotal = []
        BTotal_lowError = []

        for target in data:
            B_V = target[2] - target[4]
            if B_V > self.app.B_VMin and B_V < self.app.B_VMax:
                BTotal.append(target)
        for target in data_lowError:
            B_V = target[2] - target[4]
            if B_V > self.app.B_VMin and B_V < self.app.B_VMax:
                BTotal_lowError.append(target)

        self.NumBTotal = len(BTotal)
        self.NumBTotal_lowError = len(BTotal_lowError)

        # Ratio of Be to Be + B (total)
        self.Be_ratio = self.NumBe / self.NumBTotal
        self.Be_ratio_lowError = self.NumBe_lowError / self.NumBTotal_lowError

    def CompileBeLists(self):
        BeCandidates = []
        count = 1

        baseDate = Observations().ListDates(self.cluster)[0]
        for date in Observations().ListDates(self.cluster):
            data = np.loadtxt("../output/" + self.cluster + "/" + date + "/belist.dat")
            with fits.open("../photometry/" + self.cluster + "/" + date + "/B1.fits") as file:
                julian = file[0].header["JD"]
                binning = file[0].header["XBINNING"]
                ra = file[0].header["CRVAL1"]
                dec = file[0].header["CRVAL2"]

            if date == baseDate:
                for target in data:
                    be = []
                    be.append(target[0])   # x on image
                    be.append(target[1])   # y on image
                    be.append(target[0] * binning)   # x ref
                    be.append(target[1] * binning)   # y ref
                    be.append(self.cluster + "-WBBe" + str(count))   # identifier
                    be.append("gained")   # transient status
                    be.append(ra)   # placeholder ra
                    be.append(dec)   # placeholder dec
                    be.append(julian)   # julian date
                    be.append(date)   # date
                    be.append(count)

                    BeCandidates.append(be)
                    count += 1
            else:
                # Get x- and y- offsets
                filename = "../output/" + self.cluster + "/" + date + "/astrometry_offsets.txt"
                try:
                    offsets = np.loadtxt(filename)
                except IOError:
                    print("  Calculating coordinate offsets...")
                    offsets = Astrometry().GetOffset(self.cluster, date)
                    with open(filename, 'w') as F:
                        np.savetxt(F, offsets, fmt="%.3f")

                xOffset = offsets[0]
                yOffset = offsets[1]

                # Compare with other dates, see if target is already found
                for target in data:
                    be = []
                    newCandidate = True
                    for candidate in BeCandidates:
                        # if they refer to the same star
                        if abs(candidate[2] - (binning * target[0] + xOffset)) <= self.app.cooTol and \
                                abs(candidate[3] - (binning * target[1] + yOffset)) <= self.app.cooTol:
                            # Determine if it was gained on this date
                            transient = "gained"
                            prev_date_ind = Observations().ListDates(self.cluster).index(date) - 1
                            for item in BeCandidates:
                                if item[4] == candidate[4] and item[9] == Observations().ListDates(self.cluster)[prev_date_ind]:
                                    transient = "  --  "
                                    break

                            be.extend([target[0], target[1], binning * target[0] + xOffset, binning * target[1] + yOffset, candidate[4], transient, candidate[6], candidate[7], julian, date, candidate[10]])
                            BeCandidates.append(be)

                            newCandidate = False
                            break

                    if newCandidate:
                        be = []
                        be.append(target[0])   # x on image
                        be.append(target[1])   # y on image
                        be.append(binning * target[0] + xOffset)   # x ref
                        be.append(binning * target[1] + yOffset)   # y ref
                        be.append(self.cluster + "-WBBe" + str(count))   # identifier
                        be.append("gained")   # transient status
                        be.append(ra)   # placeholder ra
                        be.append(dec)   # placeholder dec
                        be.append(julian)   # julian date
                        be.append(date)   # date
                        be.append(count)

                        BeCandidates.append(be)
                        count += 1

        return BeCandidates

    def Summary(self):
        with open("../output/" + self.cluster + "/summary.txt", 'w') as F:
            F.write("================================================\n")
            F.write("                " + self.cluster + " Summary                  \n")
            F.write("================================================\n\n\n")

            for date in Observations().ListDates(self.cluster):
                F.write("                   " + date + "                   \n")
                F.write("------------------------------------------------\n")

                # Threshold information
                filename = "../output/" + self.cluster + "/" + date + "/thresholds.dat"
                thresholds = np.loadtxt(filename)
                F.write("Thresholds:\n")
                F.write("   Constant: " + "%.3f" % thresholds[0][1] + "     Linear: " + "%.3f" % thresholds[1][0] + ", " + "%.3f" % thresholds[1][1] + "\n\n")

                # Be candidate information
                self.Values(date)

                F.write("Be candidates:                " + str(self.NumBe) + "\n")
                F.write("Be candidates (low error):    " + str(self.NumBe_lowError) + "\n")
                F.write("B total:                      " + str(self.NumBTotal) + "\n")
                F.write("B total (low error):          " + str(self.NumBTotal_lowError) + "\n")
                F.write("Be ratio:                     " + "%.3f" % self.Be_ratio + "\n")
                F.write("Be ratio (low error):         " + "%.3f" % self.Be_ratio_lowError + "\n")
                F.write("\n\n")

        # Write compiled Be list
        data = self.CompileBeLists()
        data = sorted(data, key=lambda x: (x[10], x[9]))   # Sort by identifier then date

        identifier = [x[4] for x in data]
        transient = [x[5] for x in data]
        ra = [x[6] for x in data]
        dec = [x[7] for x in data]
        julian = [x[8] for x in data]

        filename = "../output/" + self.cluster + "/BeList.dat"
        with open(filename, 'w') as F:
            for i in range(0, len(data)):
                F.write(identifier[i] + "\t" + transient[i] + "\t" + "%.10f" % ra[i] + "\t" + "%.10f" % dec[i] + "\t" + "%.10f" % julian[i] + "\n")

        filename = "../output/" + self.cluster + "/test.txt"
        with open(filename, 'w') as F:
            F.writelines('\t'.join(str(j) for j in i) + '\n' for i in data)
