import numpy as np
import os.path
from observations import Observations
from astrometry_offsets import AstrometryOffset
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
        BType = []
        BType_lowError = []

        for target in data:
            B_V = target[2] - target[4]
            if B_V > self.app.B_VMin and B_V < self.app.B_VMax and \
                    target not in beData:
                BType.append(target)
        for target in data_lowError:
            B_V = target[2] - target[4]
            if B_V > self.app.B_VMin and B_V < self.app.B_VMax and \
                    target not in beData_lowError:
                BType_lowError.append(target)

        self.NumB = len(BType)
        self.NumB_lowError = len(BType_lowError)

        # Ratio of Be to Be + B (total)
        self.Be_ratio = self.NumBe / (self.NumBe + self.NumB)
        self.Be_ratio_lowError = self.NumBe_lowError / (self.NumBe_lowError + self.NumB_lowError)

    def CompileBeLists(self):
        BeCandidates = []
        count = 1

        baseDate = Observations.ListDates(self.cluster)[0]
        for date in Observations.ListDates(self.cluster):
            data = np.loadtxt("../output/" + self.cluster + "/" + date + "/belist.dat")
            with fits.open("../photometry/" + self.cluster + "/" + date + "/B1.fits") as file:
                julian = file[0].header["JD"]

            if date == baseDate:
                for target in data:
                    be = []
                    be.append(target[0])   # x ref
                    be.append(target[1])   # y ref
                    be.append(self.cluster + "-WBBe" + count)   # identifier
                    be.append("gained")   # transient status
                    # be.append(...)   # ra
                    # be.append(...)   # dec
                    be.append(julian)   # julian date
                    be.append(date)   # date

                    BeCandidates.append(be)
                    count += 1
            else:
                # Get x- and y- offsets
                filename = "../output/" + self.cluster + "/" + date + "/astrometry_offsets.txt"
                if not os.path.isfile(filename):
                    offsets = AstrometryOffset.GetOffset(self.cluster, self.date)
                    xOffset = offsets[0]
                    yOffset = offsets[1]
                    with open(filename, 'w') as F:
                        np.savetxt(F, offsets, fmt="%.3f")
                else:
                    offsets = np.loadtxt(filename)
                    xOffset = float(offsets[0])
                    yOffset = float(offsets[1])

                # Compare with other dates, see if target is already found
                for target in data:
                    be = []
                    newCandidate = True
                    for candidate in BeCandidates:
                        # if they refer to the same star
                        if candidate[0] - (target[0] + xOffset) <= self.app.cooTol and \
                                candidate[1] - (target[1] + yOffset) <= self.app.cooTol:
                            # Determine if it was gained on this date
                            for item in BeCandidates:
                                prev_date_ind = Observations.ListDates(self.cluster).index(date) - 1
                                if item[2] == candidate[2] and item[7] == Observations.ListDates(self.cluster)[prev_date_ind]:
                                    transient = "-"
                                else:
                                    transient = "gained"

                            be.extend(target[0], target[1], candidate[2], transient, candidate[4], candidate[5], julian, date)
                            BeCandidates.append(be)

                            newCandidate = False
                            break

                    if newCandidate:
                        be = []
                        be.append(target[0])   # x ref
                        be.append(target[1])   # y ref
                        be.append(self.cluster + "-WBBe" + count)   # identifier
                        be.append("gained")   # transient status
                        # be.append(...)   # ra
                        # be.append(...)   # dec
                        be.append(julian)   # julian date
                        be.append(date)   # date

                        BeCandidates.append(be)
                        count += 1

        return BeCandidates

    def Summary(self):
        with open("../output/" + self.cluster + "/summary.txt", 'w') as F:
            F.write("================================================\n")
            F.write("                " + self.cluster + " Summary                  \n")
            F.write("================================================\n\n\n")

            for date in Observations.ListDates(self.cluster):
                F.write("                   " + date + "                   \n")
                F.write("------------------------------------------------\n")

                # Threshold information
                filename = "../output/" + self.cluster + "/" + date + "/thresholds.dat"
                thresholds = np.loadtxt(filename)
                F.write("Thresholds:")
                F.write("   Constant: ", "{%.3f}" % thresholds[0][1], "     Linear: ", "{%.3f}" % thresholds[1][0], ", ", "{%.3f}" % thresholds[1][1], "\n")

                filename = "../output/" + self.cluster + "/" + date + "/thresholds_scaled.dat"
                thresholds = np.loadtxt(filename)
                F.write("Scaled thresholds:")
                F.write("   Constant: ", "{%.3f}" % thresholds[0][1], "     Linear: ", "{%.3f}" % thresholds[1][0], ", ", "{%.3f}" % thresholds[1][1], "\n\n")

                # Be candidate information
                self.Values(date)

                F.write("Be candidates:                ", self.NumBe)
                F.write("Be candidates (low error):    ", self.NumBe_lowError)
                F.write("Be ratio:                     ", "{%.3f}" % self.Be_ratio)
                F.write("Be ratio (low error):         ", "{%.3f}" % self.Be_ratio_lowError)

        # Write compiled Be list
        data = self.CompileBeLists()

        filename = "../output/" + self.cluster + "/BeList.dat"
        with open(filename) as F:
            np.savetxt(F, data)
