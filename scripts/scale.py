import numpy as np
from astropy.io import fits
from astrometry import Astrometry
from observations import Observations
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator, FormatStrFormatter


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
            with fits.open(filename) as file:
                binning = file[0].header["XBINNING"]
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
            isAlone = True
            for otherTarget in data:
                r = binning * np.sqrt((target[0] - otherTarget[0])**2 + (target[1] - otherTarget[1])**2)
                if r <= self.app.cooTol and target != otherTarget:
                    isAlone = False
                    break
            if not isAlone:
                data.remove(target)

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

    def Scale(self):
        """Scales the given set of data.

        Corresponds each target in the data set to others in the reference data set and
        averages magnitude differences between each set of targets to create an averaged
        magnitude offset for each filter.  Does not use either targets within a certain
        x-y distance of any others or targets that are considered outliers.

        """
        # Read data
        print("  Creating reference sample...")
        baseDate = Observations().ListDates(self.cluster)[0]   # Establish first date as scaling base
        baseBinning = self.Binning(baseDate)
        baseData = self.SetData(baseDate)

        print("  Creating variable sample...\n")
        binning = self.Binning(self.date)
        data = self.SetData(self.date)

        # Get x- and y- offsets
        filename = "../output/" + self.cluster + "/" + self.date + "/astrometry_offsets.txt"
        try:
            offsets = np.loadtxt(filename)
        except IOError:
            offsets = Astrometry().GetOffset(self.cluster, self.date, baseDate)
            with open(filename, 'w') as F:
                np.savetxt(F, offsets, fmt="%.3f")

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

        # Get statistics from scaling
        B_offset = 0
        B_std = 0
        V_offset = 0
        V_std = 0
        R_offset = 0
        R_std = 0
        H_offset = 0
        H_std = 0

        print("  \nGenerating magnitude difference histograms...")

        if B_diff != []:
            B_offset = np.mean(B_diff)
            B_std = np.std(B_diff)
            self.num_vs_mag_hist(B_diff, B_offset, B_std, "B")
        if V_diff != []:
            V_offset = np.mean(V_diff)
            V_std = np.std(V_diff)
            self.num_vs_mag_hist(V_diff, V_offset, V_std, "V")
        if R_diff != []:
            R_offset = np.mean(R_diff)
            R_std = np.std(R_diff)
            self.num_vs_mag_hist(R_diff, R_offset, R_std, "R")
        if H_diff != []:
            H_offset = np.mean(H_diff)
            H_std = np.std(H_diff)
            self.num_vs_mag_hist(H_diff, H_offset, H_std, "H-alpha")

        # Print scale information
        print("\n  Scaled with ", len(B_diff), " stars:")
        print("  B offset = " + "%.3f" % B_offset + " +/- " + "%.3f" % B_std)
        print("  V offset = " + "%.3f" % V_offset + " +/- " + "%.3f" % V_std)
        print("  R offset = " + "%.3f" % R_offset + " +/- " + "%.3f" % R_std)
        print("  H offset = " + "%.3f" % H_offset + " +/- " + "%.3f" % H_std)

        # Implement scale offsets
        filename = "../output/" + self.cluster + "/" + self.date + "/phot_" + self.app.phot_type + ".dat"
        orig_data = np.loadtxt(filename)

        for target in orig_data:
            target[2] += B_offset
            target[4] += V_offset
            target[6] += R_offset
            target[8] += H_offset

            target[3] = np.sqrt(target[3]**2 + B_std**2)
            target[5] = np.sqrt(target[5]**2 + V_std**2)
            target[7] = np.sqrt(target[7]**2 + R_std**2)
            target[9] = np.sqrt(target[9]**2 + H_std**2)

        # Write to file
        filename = "../output/" + self.cluster + "/" + self.date + "/phot_" + self.app.phot_type + "_scaled.dat"
        with open(filename, 'w') as F:
            np.savetxt(F, orig_data, fmt="%.3f")

    def num_vs_mag_hist(self, x, mean, std, filter):
        # plt.style.use('ggplot')

        # Plot main data
        plt.hist(x, bins=20, range=(mean - 3 * std, mean + 3 * std), color='#3f3f3f')
        # plt.title(filter + " Magnitude Scaling Differences")
        plt.xlabel("Magnitude Difference", fontsize=24)
        plt.ylabel("Number", fontsize=24)

        # plt.axes().xaxis.set_major_locator(MultipleLocator(0.1))
        # plt.axes().xaxis.set_major_formatter(FormatStrFormatter('%d'))
        # plt.axes().xaxis.set_minor_locator(MultipleLocator(0.025))

        # plt.axes().yaxis.set_major_locator(MultipleLocator(10))
        # plt.axes().yaxis.set_major_formatter(FormatStrFormatter('%d'))
        # plt.axes().yaxis.set_minor_locator(MultipleLocator(2.5))

        plt.axes().tick_params('both', length=6, width=2, which='major', top=True, right=True, labelsize=16)
        plt.axes().tick_params('both', length=4, width=1, which='minor', top=True, right=True)

        plt.axes().spines['top'].set_linewidth(2)
        plt.axes().spines['right'].set_linewidth(2)
        plt.axes().spines['bottom'].set_linewidth(2)
        plt.axes().spines['left'].set_linewidth(2)

        ymin, ymax = plt.ylim()
        plt.vlines(mean + std, 0, ymax, colors='#ff5151', linestyles='dashed', label='Standard Error')
        plt.vlines(mean - std, 0, ymax, colors='#ff5151', linestyles='dashed')

        plt.legend()
        plt.tight_layout()

        # Output
        filename = "../output/" + self.cluster + "/" + self.date + "/plots/num_vs_" + filter[0] + "mag_diffs_" + self.app.phot_type + ".png"
        plt.savefig(filename, dpi=300)

        plt.clf()
