import numpy as np
from sklearn.linear_model import LinearRegression


def ProcessBeFilter(cluster, date, app, scaled):
    print("\nExtracting Be candidates for " + cluster + " on " + date + "...\n")

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

    constant_threshold = ConstantAutoThreshold(data)
    linear_threshold = LinearAutoThreshold(data, app)

    thresholds = [constant_threshold, linear_threshold]

    # Write thresholds to file
    filename = 'output/' + cluster + '/' + date + '/thresholds_' + app.phot_type + '.dat'
    with open(filename, 'w') as F:
        np.savetxt(F, thresholds, fmt='%.3f')

    return thresholds


def ConstantAutoThreshold(data, iterate_limit=50):
    """Automatically determines a constant R-H threshold.

    Statisically calculates the R-H threshold by iteratively deciding which targets
    lie above three-sigma (3 times the standard deviation) of the mean R-H value.

    Args:
            iterate_limit: Maximum number of iterations that is allowed to calculate the limit.

    Returns:
            The newly calculated threshold value.

    """
    print("  Calculating constant threshold...")

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

    print("    Constant threshold found to be:")
    print("    ", threshold)
    return [0, threshold]


def LinearAutoThreshold(data, app, iterate_limit=50):
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

    points = np.column_stack((data[:, 2] - data[:, 4], data[:, 6] - data[:, 8]))   # [ [b_v[i], r_h[i]], ... ]
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

    correlation = Line(points[:, 0], points[:, 1])   # [slope, intercept, std]
    threshold = [correlation[0], correlation[1] + 3 * correlation[2]]   # [slope, intercept]

    count = 0
    while count < iterate_limit:
        new_points = []
        for point in points:
            if point[1] <= (correlation[0] * point[0] + correlation[1] + 3 * correlation[2]) and \
               point[1] >= (correlation[0] * point[0] + correlation[1] - 3 * correlation[2]):
                new_points.append(point)

        new_points = np.array(new_points)
        new_correlation = Line(new_points[:, 0], new_points[:, 1])
        threshold = [new_correlation[0], new_correlation[1] + 3 * new_correlation[2]]

        if correlation != new_correlation:
            correlation = new_correlation
        else:
            break

        count += 1

    print("    Linear threshold found to be:")
    print("    R-H = ", threshold[0], " B-V + ", threshold[1])
    return threshold


def BeFilter(app, data, thresholds):
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
