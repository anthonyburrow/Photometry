import numpy as np
from astropy.io import fits
from astrometry import Astrometry
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
from os import remove
import os.path


class Scale:
    """Creates a Scale object, used implement offsets between observations.

    Calculates average offsets between different nights and adds these offsets
    to original data to create "scaled" photometry used for comparison.

    Attributes:
            cluster: Cluster for which data to read from.
            date: Data for which data to read from.
            app: GUI application that specifies parameters.
    """

    def __init__(self, cluster, date, app):
        self.cluster = cluster
        self.date = date
        self.app = app

    def Binning(self, date):
        """Determines the bin size used during the observation.

        Reads the B1.fits image for the specified date of observation and extracts
        the binning data.  This uses x-axis binning, however x- and y- axis binning
        are typically equal.

        Args:
                date: The date of observation.

        Returns:
                The bin size as an integer.

        """
        filename = "../photometry/" + self.cluster + "/" + date + "/B1.fits"
        try:
            F = fits.getheader(filename)
            binning = F["XBINNING"]
        except IOError:
            print("\nFile does not exist:\n" + filename)
            binning = 1

        return binning

    def SetData(self, date):
        """Determines which targets are not within a given tolerance of another.

        Given an input data set, this outputs those targets which are not within a
        specified spacial tolerance in pixels.

        Args:
                date:

        Returns:
                2-dimensional array consisting of X- and Y- image coordinates, magnitudes,
                and magnitude errors.

        """
        binning = self.Binning(date)

        # Read low error data, or else read normal data
        try:
            filename = "../output/" + self.cluster + "/" + date + "/phot_" + self.app.phot_type + "_lowError.dat"
            data = np.loadtxt(filename).tolist()
        except IOError:
            try:
                filename = "../output/" + self.cluster + "/" + date + "/phot_" + self.app.phot_type + ".dat"
                data = np.loadtxt(filename).tolist()
            except IOError:
                print("\nFile does not exist:\n" + filename)
                return

        # Get rid of stars with other stars next to them
        for target in data:
            for otherTarget in data:
                r = binning * np.sqrt((target[0] - otherTarget[0])**2 + (target[1] - otherTarget[1])**2)
                if r <= self.app.cooTol and target != otherTarget:
                    data.remove(target)
                    break

        # Get rid of stars that are outliers
        try:
            filename = "../output/" + self.cluster + "/" + date + "/beList_" + self.app.phot_type + ".dat"
            filtered_data = np.loadtxt(filename).tolist()
            for target in data:
                if target in filtered_data:
                    data.remove(target)
        except IOError:
            print("  Note: Outliers were not removed from scale sample because 'beList_" + self.app.phot_type + ".dat' does not exist.")

        return data

    def Scale(self, baseDate):
        """Scales the given set of data.

        Corresponds each target in the data set to others in the reference data set and
        averages magnitude differences between each set of targets to create an averaged
        magnitude offset for each filter.  Does not use either targets within a certain
        x-y distance of any others or targets that are considered outliers.

        """
        # Reference date process
        if self.date == baseDate:
            # Set documented scale offsets to zero
            offsets = [[0, 0], [0, 0], [0, 0], [0, 0]]
            with open("../output/" + self.cluster + "/" + self.date + "/magScales.dat", 'w') as F:
                np.savetxt(F, offsets, fmt="%.3f")

            # Remove any unneeded/extra files (usually after re-scaling)
            path = "../output/" + self.cluster + "/" + self.date + "/"
            files = [
                path + "plots/num_vs_Bmag_diffs_" + self.app.phot_type + ".png",
                path + "plots/num_vs_Vmag_diffs_" + self.app.phot_type + ".png",
                path + "plots/num_vs_Rmag_diffs_" + self.app.phot_type + ".png",
                path + "plots/num_vs_Hmag_diffs_" + self.app.phot_type + ".png"
            ]
            for file in files:
                if os.path.isfile(file):
                    remove(file)

            return

        print("\nScaling data for " + self.cluster + " on " + self.date + "...\n")

        # Read data
        print("  Creating reference sample...")
        baseBinning = self.Binning(baseDate)
        baseData = self.SetData(baseDate)

        print("  Creating variable sample...")
        binning = self.Binning(self.date)
        data = self.SetData(self.date)

        print("\n")

        # Get x- and y-offsets
        offsets = Astrometry().GetOffset(self.cluster, self.date, baseDate)
        xOffset = offsets[0]
        yOffset = offsets[1]

        # Find corresponding targets
        B_diff = []
        V_diff = []
        R_diff = []
        H_diff = []

        for target in data:
            for baseTarget in baseData:
                if abs(baseBinning * baseTarget[0] - (binning * target[0] + xOffset)) <= self.app.cooTol and \
                   abs(baseBinning * baseTarget[1] - (binning * target[1] + yOffset)) <= self.app.cooTol:
                    B_diff.append(baseTarget[2] - target[2])
                    V_diff.append(baseTarget[4] - target[4])
                    R_diff.append(baseTarget[6] - target[6])
                    H_diff.append(baseTarget[8] - target[8])

        offsets = self.GetOffsets([B_diff, V_diff, R_diff, H_diff])

        # Print scale information
        print("\n  Scaled with ", len(B_diff), " stars:")
        print("  B offset = " + "%.3f" % offsets[0][0] + " +/- " + "%.3f" % offsets[0][1])
        print("  V offset = " + "%.3f" % offsets[1][0] + " +/- " + "%.3f" % offsets[1][1])
        print("  R offset = " + "%.3f" % offsets[2][0] + " +/- " + "%.3f" % offsets[2][1])
        print("  H offset = " + "%.3f" % offsets[3][0] + " +/- " + "%.3f" % offsets[3][1])

        # Output to file
        with open("../output/" + self.cluster + "/" + self.date + "/magScales.dat", 'w') as F:
            np.savetxt(F, offsets, fmt="%.3f")

    def GetOffsets(self, diffs):
        offsets = []
        for diff in diffs:
            if diff != []:
                offset = []
                offset.append(np.mean(diff))
                offset.append(np.std(diff))

                # Recalculate using only the mag differences within 3-sigma
                for target in diff:
                    if not offset[0] - 3 * offset[1] < target < offset[0] + 3 * offset[1]:
                        diff.remove(target)
                offset[0] = np.mean(diff)
                offset[1] = np.std(diff)
            else:
                offset = [0, 0]

            offsets.append(offset)

            titles = ["B", "V", "R", "H-alpha"]
            # self.num_vs_mag_hist(diff, offset[0], offset[1], titles[diffs.index(diff)])

        return offsets

    def num_vs_mag_hist(self, x, mean, std, filter):
        plt.figure(figsize=(12, 9))

        # Plot main data
        plt.hist(x, bins=20, range=(mean - 3 * std, mean + 3 * std), color='#3f3f3f')

        # plt.title(filter + " Magnitude Scaling Differences")
        plt.xlabel("Magnitude Difference", fontsize=36)
        plt.ylabel("Frequency", fontsize=36)

        plt.axes().xaxis.set_major_locator(MultipleLocator(0.1))
        plt.axes().xaxis.set_major_formatter(FormatStrFormatter('%.2f'))
        plt.axes().xaxis.set_minor_locator(MultipleLocator(0.025))

        # plt.axes().yaxis.set_major_locator(MultipleLocator(10))
        # plt.axes().yaxis.set_major_formatter(FormatStrFormatter('%d'))
        # plt.axes().yaxis.set_minor_locator(MultipleLocator(2.5))

        plt.axes().tick_params('both', length=12, width=4, which='major', top=True, right=True, direction='in', pad=6, labelsize=30)
        plt.axes().tick_params('both', length=8, width=3, which='minor', top=True, right=True, direction='in')

        plt.axes().spines['top'].set_linewidth(4)
        plt.axes().spines['right'].set_linewidth(4)
        plt.axes().spines['bottom'].set_linewidth(4)
        plt.axes().spines['left'].set_linewidth(4)

        ymin, ymax = plt.ylim()
        plt.vlines(mean + std, 0, ymax, colors='#ff5151', linestyles='dashed', label='Standard Error')
        plt.vlines(mean - std, 0, ymax, colors='#ff5151', linestyles='dashed')

        # plt.legend(fontsize=28)
        plt.tight_layout()

        # Output
        filename = "../output/" + self.cluster + "/" + self.date + "/plots/num_vs_" + filter[0] + "mag_diffs_" + self.app.phot_type + ".png"
        plt.savefig(filename)

        plt.clf()
