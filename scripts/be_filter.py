import os.path
import numpy as np
from sklearn.linear_model import LinearRegression


class BeFilter:
    """Creates a BeFilter object, holding information on determining Be candidates.

    Be candidates are extracted and output into a single respective file.  The threshold
    that determines the candidates can be calculated automatically or specified manually.

    Attributes:
            cluster:
            date:
            phot_type: Determines whether "psf" or "aperture" photometry is desired.
            auto_threshold: When true, the R-H limit will be calculated automatically.
            R_H_threshold: Manually issue a specific threshold for determining candidates.
            B_V_limits: Specifies the color range desired for the filter process.
            scaled:
    """

    def __init__(self, cluster, date, phot_type="psf", auto_threshold=True, auto_threshold_type="Linear", R_H_threshold=-3.75, B_V_limits=[-0.25, 1], scaled=False):
        self.cluster = cluster
        self.date = date
        self.phot_type = phot_type
        self.auto_threshold = auto_threshold
        self.auto_threshold_type = auto_threshold_type
        self.R_H_threshold = R_H_threshold
        self.B_V_limits = B_V_limits
        self.scaled = scaled

    def Full(self):
        """Extracts Be candidate data for every target.

        Issues the command to filter the entirety of the finalized photometry data.
        """
        with np.loadtxt("../output/" + self.cluster + "/" + self.date + "/phot_" + self.phot_type + ".dat") as data:
            if not self.scaled:
                self.Filter(data, "beList.dat")
            else:
                self.Filter(data, "beList_scaled.dat")

    def LowError(self):
        """Extracts Be candidate data for targets with lower error.

        Issues the command to filter the finalized photometry data which exhibits
        constrained error.
        """
        with np.loadtxt("../output/" + self.cluster + "/" + self.date + "/phot_" + self.phot_type + "_lowError.dat") as data:
            if not self.scaled:
                self.Filter(data, "beList_lowError.dat")
            else:
                self.Filter(data, "beList_lowError_scaled.dat")

    def Filter(self, data, output):
        """Determines which targets lie outside the threshold.

        Determines which targets lie outside the threshold and writes to a corresponding
        output.

        Args:
                data (array): Data set that is to be filtered.
                output (string): The filename for the output data.

        Returns:
                2-dimensional array consisting of X- and Y- image coordinates, magnitudes,
                and magnitude errors for each target that is filtered.

        """
        B_V = data[:, 2] - data[:, 4]
        R_H = data[:, 6] - data[:, 8]

        filtered_data = []

        if self.auto_threshold:
            if self.auto_threshold_type == "Linear":
                self.ConstantAutoThreshold()
                self.R_H_threshold = self.LinearAutoThreshold()

                for i in range(0, len(data)):
                    if R_H[i] >= (self.R_H_threshold[0] * B_V[i] + self.R_H_threshold[1]) and \
                       B_V[i] >= self.B_V_limits[0] and B_V[i] <= self.B_V_limits[1]:
                        filtered_data.append(data[i])

            elif self.auto_threshold_type == "Constant":
                self.R_H_threshold = self.ConstantAutoThreshold()
                self.LinearAutoThreshold()

                for i in range(0, len(data)):
                    if R_H[i] >= self.R_H_threshold and \
                       B_V[i] >= self.B_V_limits[0] and B_V[i] <= self.B_V_limits[1]:
                        filtered_data.append(data[i])

        # Output to file
        with open("../output/" + self.cluster + "/" + self.date + "/" + output, 'w') as F:
            for item in filtered_data:
                F.write(" ".join(item) + "\n")

        return filtered_data

    def ConstantAutoThreshold(self, iterate_limit=10):
        """Automatically determines a constant R-H threshold.

        Statisically calculates the R-H threshold by iteratively deciding which targets
        lie above three-sigma (3 times the standard deviation) of the mean R-H value.

        Args:
                iterate_limit: Maximum number of iterations that is allowed to calculate the limit.

        Returns:
                The newly calculated threshold value.

        """
        print(" Calculating constant threshold...")

        with np.loadtxt("../output/" + self.cluster + "/" + self.date + "/phot_" + self.phot_type + "_lowError.dat") as data:
            r_h = data[:, 6] - data[:, 8]

        mean = np.mean(r_h)
        std = np.std(r_h)
        threshold = mean + 3 * std

        count = 0
        while count < iterate_limit:
            new_r_h = []
            for value in r_h:
                if value <= (mean + 3 * std) and value >= (mean - 3 * std):
                    new_r_h.append(value)

            new_mean = np.mean(new_r_h)
            new_std = np.std(new_r_h)
            threshold = new_mean + 3 * new_std

            if mean != new_mean and std != new_std:
                mean = new_mean
                std = new_std
            else:
                break

            count += 1

        # Write to file

        if self.scaled:
            outputFile = "thresholds_scaled.dat"
        else:
            outputFile = "thresholds.dat"

        if os.path.isfile("../output/" + self.cluster + "/" + self.date + "/" + outputFile):
            with open("../output/" + self.cluster + "/" + self.date + "/" + outputFile) as F:
                file = F.readlines()
                file[0] = threshold   # overwrite constant threshold only
        else:
            file = [threshold, [0, threshold]]

        with open("../output/" + self.cluster + "/" + self.date + "/" + outputFile, 'w') as F:
            for item in file:
                F.write(" ".join(item) + "\n")

        print("Constant threshold found to be " + threshold)
        return threshold

    def LinearAutoThreshold(self, iterate_limit=50):
        """Automatically determines a linear R-H threshold.

        Statisically calculates the R-H threshold by iteratively deciding which targets
        lie above three-sigma (3 times the standard deviation) of a linear regression that fits data
        between the B-V limits.

        Args:
                iterate_limit: Maximum number of iterations that is allowed to calculate the limit.

        Returns:
                The newly calculated threshold line as an array consisting of slope and intercept.

        """
        print(" Calculating linear threshold...")

        with np.loadtxt("../output/" + self.cluster + "/" + self.date + "/phot_" + self.phot_type + "_lowError.dat") as data:
            points = np.column_stack((data[:, 2] - data[:, 4], data[:, 6] - data[:, 8]))   # [ [b_v[i], r_h[i]], ... ]
            for i in range(0, len(points)):
                if points[i][0] <= self.B_V_limits[0] and points[i][0] >= self.B_V_limits[1]:
                    points = np.delete(points, i)

        correlation = self.Line(points[:, 0], points[:, 1])   # [slope, intercept, std]
        threshold = [correlation[0], correlation[1] + 3 * correlation[2]]   # [slope, intercept]

        count = 0
        while count < iterate_limit:
            new_points = []
            for point in points:
                if point[1] <= (correlation[0] * point[0] + correlation[1] + 3 * correlation[2]) and \
                   point[1] >= (correlation[0] * point[0] + correlation[1] - 3 * correlation[2]):
                    new_points.append(point)

            new_correlation = self.Line(new_points[:, 0], new_points[:, 1])
            threshold = [new_correlation[0], new_correlation[1] + 3 * new_correlation[2]]

            if correlation != new_correlation:
                correlation = new_correlation
            else:
                break

            count += 1

        # Write to file

        if self.scaled:
            outputFile = "thresholds_scaled.dat"
        else:
            outputFile = "thresholds.dat"

        if os.path.isfile("../output/" + self.cluster + "/" + self.date + "/" + outputFile):
            with open("../output/" + self.cluster + "/" + self.date + "/" + outputFile) as F:
                file = F.readlines()
                file[1] = threshold   # overwrite linear threshold only
        else:
            file = [-3.5, threshold]

        with open("../output/" + self.cluster + "/" + self.date + "/" + outputFile, 'w') as F:
            for item in file:
                F.write(" ".join(item) + "\n")

        print("Linear threshold found to be: R-H = " + "{:0.2f}".format(threshold[0]) + " B-V + " + "{:0.2f}".format(threshold[1]))
        return threshold

    def Line(x, y):
        """Calcluates linear regression and standard deviation of data.

        Statisically determines the linear fit of any data, as well as the standard deviation from
        the line.

        Args:
                x: Array of x-axis values.
                y: Array of y-axis values.

        Returns:
                The regression information of the data as an array consisting of slope, intercept, and standard deviation.

        """
        line = LinearRegression.fit(x, y)

        b1 = line.coef_
        b0 = line.intercept_

        R = []
        N = len(x)
        for i in range(0, N):
            R.append((y[i] - (b1 * x[i] + b0))**2)
        s = sum(R)
        std = np.sqrt(s / (N - 1))

        correlation = [b1, b0, std]

        return correlation
