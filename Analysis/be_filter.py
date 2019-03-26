import numpy as np
from sklearn.linear_model import LinearRegression

from .low_error import ProcessLowError
from .spectral_type import AppMag
from .io import WriteClusterProperty
from .gaia import GetRParams


def ProcessBeFilter(cluster, date, app, scaled):
    """Controls full process of calculating Be candidates for a single date.

    This program calculates R-H thresholds for the specified date and uses
    these values to determine initial Be candidates for the date. Thresholds
    can be linear or constant.

    Args:
        cluster (str): Cluster to be processed.
        date (str): Date to be processed.
        app (Application): The GUI application object that controls processing.
        scaled (bool): Determines whether data is already scaled. True if scaled
                       data is to be used, and False if not.

    """
    thresholds = GetAutoThresholds(cluster, date, app, scaled)

    # Filter complete photometry file
    ext = ''
    if scaled:
        ext = '_scaled'
    filename = 'output/%s/%s/phot%s.dat' % (cluster, date, ext)

    data = np.loadtxt(filename)
    filtered_data = BeFilter(cluster, app, data, thresholds)

    # Write Be candidates file
    ext = ''
    if scaled:
        ext = '_scaled'
    filename = 'output/%s/%s/belist%s.dat' % (cluster, date, ext)

    with open(filename, 'w') as F:
        np.savetxt(F, filtered_data, fmt='%.3f')


def GetAutoThresholds(cluster, date, app, scaled):
    """Controls full process of calculating Be candidates for a single date.

    This program calculates R-H thresholds for the specified date and uses
    these values to determine initial Be candidates for the date. Thresholds
    can be linear or constant.

    Args:
        cluster (str): Cluster from which to receive thresholds.
        date (str): Date from which to receive thresholds.
        app (Application): The GUI application object that controls processing.
        scaled (bool): Determines whether data is already scaled. True if scaled
                       data is to be used, and False if not.

    Returns:
        thresholds (2-tuple): Tuple containing the constant and linear threshold
                              line objects.

    """
    # Use low-error photometry for threshold calculations
    ext = ''
    if scaled:
        ext = '_scaled'
    filename = 'output/%s/%s/phot%s.dat' % (cluster, date, ext)

    data = np.loadtxt(filename, ndmin=2)
    data = ProcessLowError(cluster, date, data)

    print("  Calculating constant threshold...")
    constant_threshold = ConstantAutoThreshold(data)

    print("\n  Calculating linear threshold...")
    linear_threshold = LinearAutoThreshold(cluster, data, app)

    print()

    thresholds = (constant_threshold, linear_threshold)

    # Write thresholds to file
    output = [[thresholds[0].slope, thresholds[0].threshold],
              [thresholds[1].slope, thresholds[1].threshold]]
    filename = 'output/%s/%s/thresholds.dat' % (cluster, date)
    with open(filename, 'w') as F:
        np.savetxt(F, output, fmt='%.3f')

    return thresholds


class ThresholdLine:
    """Object that holds information for linear model fits against given data.

    Can be used to fit data with a constant line or a linear model, and
    determines the standard error of the data around the model.

    """

    def __init__(self, slope=0, mean=0, std=0, threshold=0):
        """Instantiates a custom ThresholdLine object.

        Args:
            slope (float): Slope of the line. Defaults to 0.
            mean (float): Y-intercept of the line fit to the data. Defaults to 0.
            std (float): Standard error of data around the line. Defaults to 0.
            threshold (float): Y-intercept of the line 3-sigma above the fit
                               line. Defaults to 0.

        """
        self.slope = slope
        self.mean = mean
        self.std = std
        self.threshold = threshold

    def ConstantFit(self, y_arr):
        """Fits data to a constant line.

        Args:
            y_arr (list): Data to be fit.

        Returns:
            ThresholdLine: Returns instantiation with fit parameters.

        """
        self.mean = np.mean(y_arr)
        self.std = np.std(y_arr)

        self.threshold = self.mean + 3 * self.std

        return self

    def LinearFit(self, x_arr, y_arr):
        """Fits data to a linear model.

        Args:
            x_arr (list): X-coordinates to be fit.
            y_arr (list): Y-coordinates to be fit.

        Returns:
            ThresholdLine: Returns instantiation with fit parameters.

        """
        reg = LinearRegression()
        x_arr = x_arr.reshape(-1, 1)
        line = reg.fit(x_arr, y_arr)

        B1 = line.coef_
        B0 = line.intercept_

        R = []
        N = len(x_arr)
        for i in range(N):
            R.append((y_arr[i] - (B1 * x_arr[i] + B0))**2)
        s = sum(R)
        std = np.sqrt(s / (N - 1))

        self.mean = B0
        self.slope = B1
        self.std = std
        self.threshold = B0 + 3 * std

        return self

    def ValueFromMean(self, x):
        """Gets value from 'mean' line (line fit to data) at point x.

        Args:
            x (float): Independent (x) variable corresponding to value.

        Returns:
            float: Value on the 'mean' line at x.

        """
        return self.slope * x + self.mean

    def ValueFromThreshold(self, x):
        """Gets value from 'threshold' line (the line 3-sigma above the fit) at
           point x.

        Args:
            x (float): Independent (x) variable corresponding to value.

        Returns:
            float: Value on the 'threshold' line at x.

        """
        return self.slope * x + self.threshold


def ConstantAutoThreshold(data, iterate_limit=50):
    """Automatically calculates a constant R-H threshold.

    Statisically calculates the R-H threshold by iteratively deciding which
    targets lie above three-sigma (3 times the standard deviation) of the mean
    R-H value.

    Args:
        data (list): List of data from which to calculate the threshold.
        iterate_limit (int): Maximum number of iterations that is allowed to
                             calculate the limit.

    Returns:
        ThresholdLine: Object containing threshold line information.

    """
    r_h = data[:, 6] - data[:, 8]
    threshold = ThresholdLine().ConstantFit(r_h)

    count = 0
    while count < iterate_limit:
        new_r_h = []
        for value in r_h:
            if abs(value - threshold.mean) <= 3 * threshold.std:
                new_r_h.append(value)

        new_threshold = ThresholdLine().ConstantFit(new_r_h)

        if threshold.slope != new_threshold.slope or \
           threshold.mean != new_threshold.mean:
            threshold = new_threshold
        else:
            break

        count += 1

    print("    R-H = %.3f" % threshold.threshold)
    return threshold


def LinearAutoThreshold(cluster, data, app, iterate_limit=50):
    """Automatically determines a linear R-H threshold.

    Statisically calculates the R-H threshold by iteratively deciding which
    targets lie above three-sigma (3 times the standard deviation) of a linear
    regression that fits data between the B-V limits.

    Args:
        data (list): List of data from which to calculate the threshold.
        app (Application): The GUI application object that controls processing.
        iterate_limit (int): Maximum number of iterations that is allowed to
                             calculate the limit.

    Returns:
        ThresholdLine: Object containing threshold line information.

    """
    points = np.column_stack((data[:, 2] - data[:, 4], data[:, 6] - data[:, 8]))

    # Cut out B-V range (NOT FOR NGC7419)
    if cluster != 'NGC7419':
        n = points.shape[0]
        for i in reversed(range(n)):
            if not app.B_VMin < points[i][0] < app.B_VMax:
                points = np.delete(points, i, 0)

    # Calculate threshold
    threshold = ThresholdLine().LinearFit(points[:, 0], points[:, 1])

    count = 0
    while count < iterate_limit:
        new_points = []
        for point in points:
            if abs(point[1] - threshold.ValueFromMean(point[0])) <= \
               3 * threshold.std:
                new_points.append(point)

        new_points = np.array(new_points)
        new_threshold = ThresholdLine().LinearFit(new_points[:, 0],
                                                  new_points[:, 1])

        if threshold.slope != new_threshold.slope or \
           threshold.mean != new_threshold.mean:
            threshold = new_threshold
        else:
            break

        count += 1

    print("    R-H = %.3f B-V + %.3f" % (threshold.slope, threshold.threshold))
    return threshold


def GetVMagLimit(cluster):
    params = GetRParams(cluster)

    M_V = 0.95   # A0 limit
    m_V = AppMag(M_V, distance=params[0])

    WriteClusterProperty(cluster, 'Vmag_lim', '%.3f' % m_V)

    return m_V


def BeFilter(cluster, app, data, thresholds):
    """Determines which targets lie outside the threshold.

    Args:
        app (Application): The GUI application object that controls processing.
        data (list): Data set that is to be filtered.
        thresholds (2-tuple): Tuple containing the constant and linear threshold
                              line objects.

    Returns:
        list: Data set with removed points outside range governed by threshold
              values.

    """
    filtered_data = []

    if app.autoThresholdCheck.isChecked():
        # Automatic threshold
        if app.threshold_type == 'Linear':
            R_H_threshold = thresholds[1]
        elif app.threshold_type == 'Constant':
            R_H_threshold = thresholds[0]
    else:
        # Manual threshold
        R_H_threshold = ThresholdLine(threshold=app.threshold)

    Vlim = GetVMagLimit(cluster)

    for target in data:
        b_v = target[2] - target[4]
        r_h = target[6] - target[8]
        r_herr = np.sqrt(target[7]**2 + target[9]**2)

        if r_h - 3 * r_herr >= R_H_threshold.ValueFromThreshold(b_v) and \
           app.B_VMin <= b_v <= app.B_VMax and \
           target[4] < Vlim:
            filtered_data.append(target)

    return filtered_data
