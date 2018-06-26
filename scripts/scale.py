import os.path
import numpy as np
from astropy.io import fits
from astrometry_offset import AstrometryOffset
from observations import Observations


class Scale:
    """Creates a Scale object, used implement offsets between observations.

    Calculates average offsets between different nights and adds these offsets
    to original data to create "scaled" photometry used for comparison.

    Attributes:
            cluster: Cluster whose data is to be scaled.
            phot_type: Determines whether "psf" or "aperture" photometry is desired.
            coo_tol: Spacial tolerance for matching stars between observations.
    """

    def __init__(self, cluster, date, app):
        self.cluster = cluster
        self.date = date
        self.app = app

        self.baseBin = 1
        self.baseData = self.BaseData()

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
        if os.path.isfile(filename):
            with fits.open(filename) as file:
                Bin = file[0].header["XBINNING"]

        return Bin

    def Filter(data, filterTol=20):
        """Determines which targets are not within a given tolerance of another.

        Given an input data set, this outputs those targets which are not within a
        specified spacial tolerance in pixels.

        Args:
                data: 2-dimensional array of data to be filtered.
                filterTol: Spacial tolerance to be used, given in pixels.

        Returns:
                2-dimensional array consisting of X- and Y- image coordinates, magnitudes,
                and magnitude errors.

        """
        filtered_data = []

        for target in data:
            isAlone = True
            for otherTarget in [x for x in data if x != target]:
                r = np.sqrt((target[0] - otherTarget[0])**2 + (target[1] - otherTarget[1])**2)
                if r > filterTol:
                    isAlone = False
                    break
            if isAlone:
                filtered_data.append(target)

        return filtered_data

    def BaseData(self):
        """Provides the reference data set used for the cluster.

        Selects the first night of observation as a reference data set to which all other
        observations for a cluster are scaled.

        Returns:
                2-dimensional array consisting of X- and Y- image coordinates, magnitudes,
                and magnitude errors for the reference data set.

        """
        date = Observations.ListDates(self.cluster)[0]  # Establish first date as scaling base

        filename = "../ouput/" + self.cluster + "/" + date + "/phot_" + self.app.phot_type + "_scaled.dat"
        if not os.path.isfile(filename):
            # Create base data and Filter (spacial filter)

            data = np.loadtxt("../output/" + self.cluster + "/" + date + "/phot_" + self.app.phot_type + ".dat")
            data = self.Filter(data, 20)
            self.baseBin = self.Binning(date)

            # Remove outliers
            with open("../output/" + self.cluster + "/" + date + "/beList.dat") as filtered_data:
                for target in data:
                    if target in filtered_data:
                        data.remove(target)

            # Write to file
            with open("../ouput/" + self.cluster + "/" + date + "/phot_" + self.app.phot_type + "_scaled.dat", 'w') as F:
                for item in data:
                    F.write(" ".join(item) + "\n")
        else:
            data = np.loadtxt(filename)

        return data

    def Scale(self):
        """Scales the given set of data.

        Corresponds each target in the data set to others in the reference data set and
        averages magnitude differences between each set of targets to create an averaged
        magnitude offset for each filter.  Does not use either targets within a certain
        x-y distance of any others or targets that are considered outliers.

        Args:
                date: Observation date of the data to be compared to reference data.

        Returns:
                Original 2-dimensional array of data with offset correction implemented.

        """
        # Load data and Filter (spacial filter)
        orig_data = np.loadtxt("../output/" + self.cluster + "/" + self.date + "/phot_" + self.app.phot_type + ".dat")
        data = self.Filter(orig_data, 20)
        binning = self.Binning(self.date)

        # Remove outliers
        with open("../output/" + self.cluster + "/" + self.date + "/beList.dat") as filtered_data:
            for target in data:
                if target in filtered_data:
                    data.remove(target)

        # Get x- and y- offsets
        filename = "../output/" + self.cluster + "/" + self.date + "/astrometry_offsets.txt"
        if not os.path.isfile(filename):
            offsets = AstrometryOffset.GetOffset(self.cluster, self.date)
            xOffset = offsets[0]
            yOffset = offsets[1]
            with open(filename, 'w') as F:
                F.write(xOffset + "\n" + yOffset)
        else:
            with open(filename) as F:
                offsets = F.readlines()
                xOffset = float(offsets[0])
                yOffset = float(offsets[1])

        # Find corresponding targets
        B_diff = []
        V_diff = []
        R_diff = []
        H_diff = []

        for target in data:
            for baseTarget in self.baseData:
                if abs(self.baseBin * baseTarget[0] - (binning * target[0] + xOffset)) <= self.app.cooTol and \
                   abs(self.baseBin * baseTarget[1] - (binning * target[1] + yOffset)) <= self.app.cooTol:
                    B_diff.append(baseTarget[2] - target[2])
                    V_diff.append(baseTarget[4] - target[4])
                    R_diff.append(baseTarget[6] - target[6])
                    H_diff.append(baseTarget[8] - target[8])

        # Get statistics from scaling
        B_offset = np.mean(B_diff)
        V_offset = np.mean(V_diff)
        R_offset = np.mean(R_diff)
        H_offset = np.mean(H_diff)

        B_std = np.std(B_diff)
        V_std = np.std(V_diff)
        R_std = np.std(R_diff)
        H_std = np.std(H_diff)

        # Write scale information
        filename = "../output/" + self.cluster + "/scale_output.dat"
        with open(filename, "a") as F:
            F.write(self.cluster + "\n\n")
            F.write("B offset = " + str(B_offset) + " +/- " + str(B_std) + "\n")
            F.write("V offset = " + str(V_offset) + " +/- " + str(V_std) + "\n")
            F.write("R offset = " + str(R_offset) + " +/- " + str(R_std) + "\n")
            F.write("H offset = " + str(H_offset) + " +/- " + str(H_std) + "\n")
            F.write("\n")

        # Implement scale offsets
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
        filename = "../ouput/" + self.cluster + "/" + self.date + "/phot_" + self.app.phot_type + "_scaled.dat"
        with open(filename, 'w') as F:
            for item in orig_data:
                F.write(" ".join(item) + "\n")
