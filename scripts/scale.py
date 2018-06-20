import numpy as np
from astropy.io import fits
from astrometry_offset import AstrometryOffset


class Scale:
    """Creates a Scale object, used implement offsets between observations.

    Calculates average offsets between different nights and adds these offsets
    to original data to create "scaled" photometry used for comparison.

    Attributes:
            cluster: Cluster whose data is to be scaled.
            phot_type: Determines whether "psf" or "aperture" photometry is desired.
            coo_tol: Spacial tolerance for matching stars between observations.
    """

    def __init__(self, cluster, phot_type="psf", coo_tol=10):
        self.cluster = cluster
        self.phot_type = phot_type
        self.coo_tol = coo_tol

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
        with fits.open("../photometry/" + self.cluster + "/" + date + "/B1.fits") as file:
            Bin = file[0].header["XBINNING"]  # Need to check

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
        # Create base data and Filter (spacial filter)
        with open("../photometry/" + self.cluster + "/obs_dates.txt") as F:  # Establish first date as scaling base
            date = F.readline()
        data = np.loadtxt("../output/" + self.cluster + "/" + date + "/phot_" + self.phot_type + ".dat")
        data = self.Filter(data, 20)
        self.baseBin = self.Binning(date)

        # Remove outliers
        with open("../output/" + self.cluster + "/" + date + "/beList.dat") as filtered_data:
            for target in data:
                if target in filtered_data:
                    data.remove(target)

        return data

    def Scale(self, date):  # TODO: Call Scale and automate offsets
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
        # Create base data and Filter (spacial filter)
        print("  Scaling all data for " + self.cluster)

        orig_data = np.loadtxt("../output/" + self.cluster + "/" + date + "/phot_" + self.phot_type + ".dat")
        data = self.Filter(orig_data, 20)
        binning = self.Binning(date)

        # Remove outliers
        with open("../output/" + self.cluster + "/" + date + "/beList.dat") as filtered_data:
            for target in data:
                if target in filtered_data:
                    data.remove(target)

        B_diff = []
        V_diff = []
        R_diff = []
        H_diff = []

        offsets = AstrometryOffset.GetOffset(self.cluster, date)
        xOffset = offsets[0]
        yOffset = offsets[1]

        for target in data:
            for baseTarget in self.baseData:
                if abs(self.baseBin * baseTarget[0] - (binning * target[0] + xOffset)) <= self.coo_tol and \
                   abs(self.baseBin * baseTarget[1] - (binning * target[1] + yOffset)) <= self.coo_tol:
                    B_diff.append(baseTarget[2] - target[2])
                    V_diff.append(baseTarget[4] - target[4])
                    R_diff.append(baseTarget[6] - target[6])
                    H_diff.append(baseTarget[8] - target[8])

        B_offset = np.mean(B_diff)
        V_offset = np.mean(V_diff)
        R_offset = np.mean(R_diff)
        H_offset = np.mean(H_diff)

        B_std = np.std(B_diff)
        V_std = np.std(V_diff)
        R_std = np.std(R_diff)
        H_std = np.std(H_diff)

        # Write scale information
        output = "../output/" + self.cluster + "/scale_output.dat"
        with open(output, "a") as F:
            F.write(self.cluster + "\n\n")
            F.write("B offset = " + str(B_offset) + " +/- " + str(B_std) + "\n")
            F.write("V offset = " + str(V_offset) + " +/- " + str(V_std) + "\n")
            F.write("R offset = " + str(R_offset) + " +/- " + str(R_std) + "\n")
            F.write("H offset = " + str(H_offset) + " +/- " + str(H_std) + "\n")
            F.write("\n")

        for target in orig_data:
            target[2] += B_offset
            target[4] += V_offset
            target[6] += R_offset
            target[8] += H_offset

            target[3] = np.sqrt(target[3]**2 + B_std**2)
            target[5] = np.sqrt(target[5]**2 + V_std**2)
            target[7] = np.sqrt(target[7]**2 + R_std**2)
            target[9] = np.sqrt(target[9]**2 + H_std**2)

        return orig_data
