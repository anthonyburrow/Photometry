import numpy as np
from sklearn.linear_model import LinearRegression


class BeFilter:
    """Creates a BeFilter object, holding information on determining Be candidates.

    Be candidates are extracted and output into a single respective file.  The threshold
    that determines the candidates can be calculated automatically or specified manually.

    Attributes:
            cluster:
            date:
            app:
            scaled:
    """

    def __init__(self, cluster, date, app, scaled):
        self.cluster = cluster
        self.date = date
        self.app = app
        self.scaled = scaled

    def Process(self):
        if self.scaled:
            filename = "../output/" + self.cluster + "/" + self.date + "/phot_" + self.app.phot_type + "_scaled.dat"
        else:
            filename = "../output/" + self.cluster + "/" + self.date + "/phot_" + self.app.phot_type + ".dat"

        try:
            data = np.loadtxt(filename)
            self.Filter(data, "beList_" + self.app.phot_type + ".dat")
        except IOError:
            print("\nFile does not exist:\n" + filename)
            return

        filename = "../output/" + self.cluster + "/" + self.date + "/phot_" + self.app.phot_type + "_lowError.dat"
        try:
            data = np.loadtxt(filename)
            self.Filter(data, "beList_" + self.app.phot_type + "_lowError.dat")
        except IOError:
            print("\nFile does not exist:\n" + filename)
            return

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

        # Automatic threshold
        if self.app.autoThresholdCheck.isChecked():
            if self.app.threshold_type == "Linear":
                self.ConstantAutoThreshold()
                R_H_threshold = self.LinearAutoThreshold()

                for i in range(0, len(data)):
                    if R_H[i] >= (R_H_threshold[0] * B_V[i] + R_H_threshold[1]) and \
                       B_V[i] >= self.app.B_VMin and B_V[i] <= self.app.B_VMax and \
                       data[i][4] <= 13.51:
                        filtered_data.append(data[i])

            elif self.app.threshold_type == "Constant":
                R_H_threshold = self.ConstantAutoThreshold()
                self.LinearAutoThreshold()

                for i in range(0, len(data)):
                    if R_H[i] >= R_H_threshold and \
                       B_V[i] >= self.app.B_VMin and B_V[i] <= self.app.B_VMax and \
                       data[i][4] <= 13.51:
                        filtered_data.append(data[i])
        # Manual threshold
        else:
            for i in range(0, len(data)):
                if R_H[i] >= R_H_threshold and \
                   B_V[i] >= self.app.B_VMin and B_V[i] <= self.app.B_VMax and \
                   data[i][4] <= 13.51:
                    filtered_data.append(data[i])

        # Output to file
        filename = "../output/" + self.cluster + "/" + self.date + "/" + output
        with open(filename, 'w') as F:
            np.savetxt(F, filtered_data, fmt='%.3f')

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
        print("  Calculating constant threshold...")

        filename = "../output/" + self.cluster + "/" + self.date + "/phot_" + self.app.phot_type + "_lowError.dat"
        try:
            data = np.loadtxt(filename)
            r_h = data[:, 6] - data[:, 8]
        except IOError:
            print("\nFile does not exist:\n" + filename)
            return

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
        filename = "../output/" + self.cluster + "/" + self.date + "/thresholds_" + self.app.phot_type + ".dat"
        try:
            file = np.loadtxt(filename)
            file[0] = [0, threshold]   # overwrite constant threshold only
        except IOError:
            file = np.array([[0, threshold], [0, threshold]])

        with open(filename, 'w') as F:
            np.savetxt(F, file, fmt='%.3f')

        print("    Constant threshold found to be:")
        print("    ", threshold)
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
        print("  Calculating linear threshold...")

        filename = "../output/" + self.cluster + "/" + self.date + "/phot_" + self.app.phot_type + "_lowError.dat"
        try:
            data = np.loadtxt(filename)
            points = np.column_stack((data[:, 2] - data[:, 4], data[:, 6] - data[:, 8]))   # [ [b_v[i], r_h[i]], ... ]
            for i in range(0, len(points)):
                if points[i][0] <= self.app.B_VMin and points[i][0] >= self.app.B_VMax:
                    points = np.delete(points, i)
        except IOError:
            print("\nFile does not exist:\n" + filename)
            return

        correlation = self.Line(points[:, 0], points[:, 1])   # [slope, intercept, std]
        threshold = [correlation[0], correlation[1] + 3 * correlation[2]]   # [slope, intercept]

        count = 0
        while count < iterate_limit:
            new_points = []
            for point in points:
                if point[1] <= (correlation[0] * point[0] + correlation[1] + 3 * correlation[2]) and \
                   point[1] >= (correlation[0] * point[0] + correlation[1] - 3 * correlation[2]):
                    new_points.append(point)

            new_points = np.array(new_points)
            new_correlation = self.Line(new_points[:, 0], new_points[:, 1])
            threshold = [new_correlation[0], new_correlation[1] + 3 * new_correlation[2]]

            if correlation != new_correlation:
                correlation = new_correlation
            else:
                break

            count += 1

        # Write to file
        filename = "../output/" + self.cluster + "/" + self.date + "/thresholds_" + self.app.phot_type + ".dat"
        try:
            file = np.loadtxt(filename)
            file[1] = threshold   # overwrite linear threshold only
        except IOError:
            file = np.array([[0, -3.5], threshold])

        with open(filename, 'w') as F:
            np.savetxt(F, file, fmt='%.3f')

        print("    Linear threshold found to be:")
        print("    R-H = ", threshold[0], " B-V + ", threshold[1])
        return threshold

    def Line(self, x, y):
        """Calcluates linear regression and standard deviation of data.

        Statisically determines the linear fit of any data, as well as the standard deviation from
        the line.

        Args:
                x: Array of x-axis values.
                y: Array of y-axis values.

        Returns:
                The regression information of the data as an array consisting of slope, intercept, and standard deviation.

        """
        reg = LinearRegression()
        x = x.reshape(-1, 1)
        line = reg.fit(x, y)

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
