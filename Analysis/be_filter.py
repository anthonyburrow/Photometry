import numpy as np
from sklearn.linear_model import LinearRegression


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
    filename = 'output/' + cluster + '/' + date + '/phot'
    if scaled:
        filename += '_scaled'
    filename += '_' + app.phot_type + '.dat'

    try:
        data = np.loadtxt(filename)
        filtered_data = BeFilter(app, data, thresholds)
    except IOError:
        print("\nFile does not exist:\n" + filename)
        return

    # Write Be candidates file
    filename = 'output/' + cluster + '/' + date + '/beList'
    if scaled:
        filename += '_scaled'
    filename += '_' + app.phot_type + '.dat'

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
        list: List containing constant threshold and linear threshold
              information.

    """
    # Use low-error photometry for threshold calculations
    filename = 'output/' + cluster + '/' + date + '/phot'
    if scaled:
        filename += '_scaled'
    filename += '_' + app.phot_type + '_lowError.dat'

    try:
        data = np.loadtxt(filename)
    except IOError:
        print("\nFile does not exist:\n" + filename)
        return

    print("  Calculating constant threshold...")
    constant_threshold = ConstantAutoThreshold(data)

    print("\n  Calculating linear threshold...")
    linear_threshold = LinearAutoThreshold(data, app)

    print()

    thresholds = [constant_threshold, linear_threshold]

    # Write thresholds to file
    filename = 'output/' + cluster + '/' + date + \
               '/thresholds_' + app.phot_type + '.dat'
    with open(filename, 'w') as F:
        np.savetxt(F, thresholds, fmt='%.3f')

    return thresholds


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
        list: List containing a constant zero slope and the calculated
              y-intercept representing the constant threshold.

    """
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

    print("    R-H = %.3f" % threshold)
    return [0, threshold]


def LinearAutoThreshold(data, app, iterate_limit=50):
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
        list: List containing a the calculated slope and y-intercept of the
              fit line.

    """
    points = np.column_stack((data[:, 2] - data[:, 4], data[:, 6] - data[:, 8]))

    n = len(points.tolist())
    for i in reversed(range(0, n)):
        if points[i][0] <= app.B_VMin or points[i][0] >= app.B_VMax:
            points = np.delete(points, i, 0)

    def Line(x, y):
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

    # [slope, intercept, std]
    correlation = Line(points[:, 0], points[:, 1])
    # [slope, intercept]
    threshold = [correlation[0], correlation[1] + 3 * correlation[2]]

    count = 0
    while count < iterate_limit:
        new_points = []
        for point in points:
            if point[1] <= (correlation[0] * point[0] +
                            correlation[1] + 3 * correlation[2]) and \
               point[1] >= (correlation[0] * point[0] +
                            correlation[1] - 3 * correlation[2]):
                new_points.append(point)

        new_points = np.array(new_points)
        new_correlation = Line(new_points[:, 0], new_points[:, 1])
        threshold = [new_correlation[0], new_correlation[1] +
                     3 * new_correlation[2]]

        if correlation != new_correlation:
            correlation = new_correlation
        else:
            break

        count += 1

    print("    R-H = %.3f" % threshold[0][0] + " B-V + %.3f" % threshold[1][0])

    return threshold


def BeFilter(app, data, thresholds):
    """Determines which targets lie outside the threshold.

    Args:
        app (Application): The GUI application object that controls processing.
        data (list): Data set that is to be filtered.
        thresholds (list): List containing constant threshold and linear
                           threshold information.

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
        R_H_threshold = [0, app.threshold]

    for target in data:
        b_v = target[2] - target[4]
        r_h = target[6] - target[8]
        if r_h >= (R_H_threshold[0] * b_v + R_H_threshold[1]) and \
           app.B_VMin <= b_v <= app.B_VMax and \
           target[4] <= 13.51:
            filtered_data.append(target)

    return filtered_data
