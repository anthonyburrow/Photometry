from mypytools.math.gauss import gauss, bimodal, hist_fit_bimodal, gauss_2d, hist_fit
import numpy as np
from scipy.optimize import curve_fit
from math import isnan

import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.ticker import MultipleLocator, FormatStrFormatter

from .io import Binning, GetWCS, WriteClusterProperty, ReadClusterProperty
from .observations import ListDates


def ProcessGaia(cluster):
    """Controls full process in calculating distance parameters.

    Distance parameters are used to determine cluster membership. Calculates
    distances parameters radially by using GAIA distance information. Calculates
    x-y pixel distances based on Gaussian fit of observed data.

    Args:
        cluster (str): The cluster from which to process distances.

    """
    print('  Fitting radial distribution...')
    ProcessRadialDistances(cluster)

    print('\n  Fitting x-y distributions...')
    ProcessXYDistances(cluster)

    print('\n  Fitting proper motion distributions...\n')
    ProcessProperMotions(cluster)

    print('  Separating targets inside and outside the cluster...\n')
    for date in ListDates(cluster):
        SeparatePhotometry(cluster, date)


def ProcessRadialDistances(cluster):
    """Controls full process in calculating radial distance parameters.

    Args:
        cluster (str): The cluster from which to process distances.

    """
    all_distances = []
    sample_params = []

    # Fit radial distance distribution for each date
    for date in ListDates(cluster):
        data = GetGAIAInfo(cluster, date)
        distances = [x[0] for x in data]
        all_distances.extend(distances)

        cluster_params = RFit(cluster, date, distances)
        sample_params.append([len(distances), cluster_params[0],
                              cluster_params[1]**2, cluster_params[2]])

    # Statistically average all dates
    RPopulationValues(cluster, all_distances, sample_params)


def ProcessXYDistances(cluster):
    """Controls full process in calculating x-y distance parameters.

    Args:
        cluster (str): The cluster from which to process distances.

    """
    for date in ListDates(cluster):
        XYFit(cluster, date)


def ProcessProperMotions(cluster):
    for date in ListDates(cluster):
        PMFit(cluster, date)


def GetGAIAInfo(cluster, date, dist_range=(0, 15)):
    """Retrieves simplified information from data exported by GAIA.

    Args:
        cluster (str): The cluster from which to retrieve information.
        date (str): The date from which to retrieve information.
        dist_range (2-tuple): Distance (in Mpc) between which information is
                              collected. Default is between 0 and 15 Mpc.

    Returns:
        list: List of the distance and RA/Dec for each target from GAIA data.

    """
    filename = 'photometry/' + cluster + '/' + date + '/phot_dists.csv'
    try:
        data = np.genfromtxt(filename, skip_header=1, usecols=(10, 98, 99, 13, 15),
                             delimiter=',')   # parallax, ra, dec, pmra, pmdec
    except IOError:
        print("Note: Data on distances not found for %s" % date)

    # Don't use negative parallax:
    data = np.array([x for x in data if x[0] > 0])
    for target in data:
        target[0] = 1 / target[0]   # convert parallax [mas] to distance [kpc]
    data = set(tuple(x) for x in data)
    data = [list(x) for x in data]   # list of unique data

    data = [x for x in data if dist_range[0] < x[0] < dist_range[1]]

    return data


def RFit(cluster, date, distances):
    """Calculates radial distance parameters.

    Parameters are calculated by histogram fits using the distance information
    obtained through GAIA. Fits are bimodal Gaussians, where the cluster
    corresponds to the peak closest to the median (with the largest amplitude).

    Args:
        cluster (str): The cluster for which fits are calculated.
        date (str): The date for which fits are calculated.
        distances (list): List of distance information from GAIA data.

    Returns:
        list: List containing the mean, standard deviation, and amplitude of
              the calculated fit for the cluster.

    """
    low_lim = 0
    up_lim = 10

    plt.style.use('researchpaper')
    fig, ax = plt.subplots()

    y, x, _ = ax.hist(distances, bins=75, range=[low_lim, up_lim],
                      color='#3f3f3f')
    x = (x[1:] + x[:-1]) / 2

    ax.set_xlabel('Distance (kpc)')
    ax.set_ylabel('Frequency')

    ax.set_xlim(low_lim, up_lim)

    median_distance = np.median(distances)
    ax.axvline(x=median_distance, linestyle='dashed', linewidth=1.5,
               color='#ff8484', label='Median = %.3f' % median_distance)

    params = hist_fit_bimodal(distances, p0=[0.8, 0.5, 3.0, 0.5], hist_bins=75,
                              hist_range=[low_lim, up_lim],
                              fit_bounds=((low_lim, 0, 0, low_lim, 0, 0),
                                          (up_lim, up_lim - low_lim, np.inf,
                                           up_lim, up_lim - low_lim, np.inf)))

    # Use curve closer to median (would use amplitude comparison instead of
    # mean but it was acting wonky)
    if abs(params[0] - median_distance) < abs(params[3] - median_distance):
        cluster_params = params[0:3]
        curve = 'G1'
    else:
        cluster_params = params[3:6]
        curve = 'G2'

    ax.axvline(x=cluster_params[0], linestyle='dashed', linewidth=1.5,
               color='#54ff79', label=r'$\mu$' + '(' + curve + ') = ' +
                                      '%.3f' % cluster_params[0])
    ax.axvline(x=cluster_params[0] + 3 * cluster_params[1], linestyle='dashed',
               linewidth=1.5, color='#5ec1ff')
    ax.axvline(x=cluster_params[0] - 3 * cluster_params[1], linestyle='dashed',
               linewidth=1.5, color='#5ec1ff')

    ax.plot(x, bimodal(x, *params), color='#ff5454', lw=4)
    ax.plot(x, gauss(x, *params[0:3]), color='#5495ff', ls='--', lw=4,
            label='G1')
    ax.plot(x, gauss(x, *params[3:6]), color='#54ff79', ls='--', lw=4,
            label='G2')

    spine_lw = 4
    [ax.spines[axis].set_linewidth(spine_lw)
     for axis in ['top', 'bottom', 'left', 'right']]

    ax.legend()

    filename = 'output/' + cluster + '/' + date + '/Rdist_distribution.png'
    fig.savefig(filename)

    plt.close('all')

    return cluster_params


def RPopulationValues(cluster, all_distances, sample_params):
    """Calculates radial distance population parameters.

    Population parameters are calculated statistically using the sample
    parameters calculated for each date. These parameters are those used
    in the final determination of cluster membership in the radial direction.

    Args:
        cluster (str): The cluster for which parameters are calculated.
        all_distances (list): List of all distances from each date.
        sample_params (list): List of all parameters from each date.

    """
    k = len(sample_params)   # Number of observations (dates)

    if k != 0:
        population_mean = np.mean(sample_params, axis=0)[1]

        z = 0
        for param in sample_params:
            z += (param[0] - 1) * param[2]
        variance = z / (np.sum(sample_params, axis=0)[0] - k)
        population_std = np.sqrt(variance)

        print("    Cluster population found to be at %.3f kpc +/- %.3f kpc" %
              (population_mean, population_std))
    else:
        print("No distance data found for %s.  Using d = 3.0 kpc +/- 1.0 kpc" %
              cluster)
        population_mean = 3.0
        population_std = 1.0

    WriteClusterProperty(cluster, 'distance', '%.3f' % population_mean)
    WriteClusterProperty(cluster, 'distance_sigma', '%.3f' % population_std)

    # Plot population data
    low_lim = 0
    up_lim = 10
    bins = 100

    plt.style.use('researchpaper')
    fig, ax = plt.subplots()

    y, x, _ = ax.hist(all_distances, bins=100, range=[low_lim, up_lim],
                      color='#3f3f3f')
    x = (x[1:] + x[:-1]) / 2  # get centers
    plt.cla()

    ax.bar(x, y / k, (up_lim - low_lim) / bins, align='center', color='#3f3f3f')

    ax.set_xlabel('Distance (kpc)')
    ax.set_ylabel('Frequency')

    ax.set_xlim(low_lim, up_lim)

    def superpose(x, params):
        g = np.zeros_like(x)
        for param in params:
            g += gauss(x, *param[1:])
        g = g / len(params)   # Normalize
        return g

    ax.axvline(x=population_mean, linestyle='dashed', linewidth=1.5,
               color='#54ff79', label=r'$\mu$' + ' = ' +
                                      '%.3f' % population_mean)
    ax.axvline(x=population_mean + 3 * population_std, linestyle='dashed',
               linewidth=1.5, color='#dc5eff', label='3-sigma')
    ax.axvline(x=population_mean - 3 * population_std, linestyle='dashed',
               linewidth=1.5, color='#dc5eff')

    # Weighted sum of gaussians
    ax.plot(x, superpose(x, sample_params), color='#f44242', ls='-', lw=4)

    spine_lw = 4
    [ax.spines[axis].set_linewidth(spine_lw)
     for axis in ['top', 'bottom', 'left', 'right']]

    ax.legend()

    filename = 'output/' + cluster + '/Rdist_pop_distribution.png'
    fig.savefig(filename)

    plt.close('all')


def XYFit(cluster, date):
    """Calculates x-y distance parameters.

    Parameters are calculated by histogram fits using the image coordinates
    of targets found for the specified date. Date-specific parameters are
    output to files.

    Args:
        cluster (str): The cluster for which fits are calculated.
        date (str): The date for which fits are calculated.

    """
    filename = 'output/' + cluster + '/' + date + '/phot.dat'
    data = np.loadtxt(filename).tolist()

    binning = Binning(cluster, date)

    x_arr = [x[0] for x in data]
    y_arr = [x[1] for x in data]

    fig = plt.figure(figsize=(16, 9))

    gs = GridSpec(2, 2)

    ax = fig.add_subplot(gs[:, 0])
    ax1 = fig.add_subplot(gs[0, 1])
    ax2 = fig.add_subplot(gs[1, 1])

    bins = 40
    h_x, x, _ = ax1.hist(x_arr, bins=bins, color='#3f3f3f')
    h_y, y, _ = ax2.hist(y_arr, bins=bins, color='#3f3f3f')

    x = np.array((x[1:] + x[:-1]) / 2)
    y = np.array((y[1:] + y[:-1]) / 2)

    def composite(x, B1, B0, mu, sigma, A):
        return gauss(x, mu, sigma, A) + B1 * x + B0

    p0_x = [0, 25, 2000 / binning, 500 / binning, max(h_x)]
    p0_y = [0, 25, 2000 / binning, 500 / binning, max(h_y)]

    params_x, cov = curve_fit(composite, x, h_x, p0=p0_x)
    params_y, cov = curve_fit(composite, y, h_y, p0=p0_y)

    ax1.plot(x, gauss(x, *params_x[-3:]) + x * params_x[0] + params_x[1],
             color='#5495ff', ls='--', lw=4, label='sigma = ' +
                                                   '%.3f' % params_x[3])
    ax2.plot(y, gauss(y, *params_y[-3:]) + y * params_y[0] + params_y[1],
             color='#5495ff', ls='--', lw=4, label='sigma = ' +
                                                   '%.3f' % params_y[3])

    ax1.set_xlim([0, 4096 / binning])
    ax1.set_ylim([0, 50])
    ax1.set_title('X', fontsize=14)

    ax2.set_xlim([0, 4096 / binning])
    ax2.set_ylim([0, 60])
    ax2.set_title('Y', fontsize=14)

    ax1.legend(fontsize=14)
    ax2.legend(fontsize=14)

    # Plot 2d gaussian
    h, x, y, _ = ax.hist2d(x_arr, y_arr, bins=bins)
    x = np.array((x[1:] + x[:-1]) / 2)
    y = np.array((y[1:] + y[:-1]) / 2)
    x, y = np.meshgrid(x, y)
    h = h.ravel()

    def composite(coo, a, b, c, mux, muy, sigma, A):
        x, y = coo
        g = gauss_2d([x, y], mux, muy, sigma, A) + a * x + b * y + c
        return g.ravel()

    p0 = [0, 0, 1, 2000 / binning, 2000 / binning, 500 / binning, max(h)]

    params, cov = curve_fit(composite, (x, y), h, p0=p0)

    ax.set_xlabel('X', fontsize=14)
    ax.set_ylabel('Y', fontsize=14)

    filename = 'output/' + cluster + '/' + date + '/XYdist_distribution.png'
    fig.savefig(filename)

    plt.close('all')

    # Output params to file
    params = params[3:5].tolist()   # mux, muy

    if params_x[3] > params_y[3]:
        params.append(params_x[3])   # sigma
    else:
        params.append(params_y[3])   # sigma

    filename = 'output/' + cluster + '/' + date + '/XYdist_params.dat'
    with open(filename, 'w') as F:
        np.savetxt(F, params, fmt='%.3f')


def GetRParams(cluster):
    """Gets radial distance parameters from output file.

    Args:
        cluster (str): The cluster for which parameters are retrieved.

    Returns:
        list: List containing cluster population mean and standard deviation.
              Returns None if output file does not exist.

    """
    params = ReadClusterProperty(cluster, ('distance', 'distance_sigma'))

    return [float(param) for param in params]


def GetXYParams(cluster, date):
    """Gets radial distance parameters from output file.

    Args:
        cluster (str): The cluster for which parameters are retrieved.
        date (str): The date for which parameters are retrieved.

    Returns:
        list: List containing cluster x- and y-coordinate means and standard
              deviation. Returns None if output file does not exist.

    """
    filename = 'output/' + cluster + '/' + date + '/XYdist_params.dat'
    try:
        params = np.loadtxt(filename)
    except IOError:
        print("'%s' does not exist." % filename)
        return

    return params.tolist()


def GetDistanceOutliers(cluster, date):
    """Gets list of RA/Decs for cluster outliers.

    Args:
        cluster (str): The cluster for which outliers are retrieved.
        date (str): The date for which outliers are retrieved.

    Returns:
        list: List containing RA/Decs for all targets that are found to be
              outside the cluster area.

    """
    # Get radial outliers
    rmu, rstd = GetRParams(cluster)
    data = GetGAIAInfo(cluster, date)

    # Get RA/Decs of stars possibly in the cluster
    radecs_in = []
    data_in = [x for x in data if abs(x[0] - rmu) < 3 * rstd]
    for target in data_in:
        radec = (target[1], target[2])
        if radec not in radecs_in:
            radecs_in.append(radec)

    # Get RA/Decs of stars with no possibility of being in cluster
    outliers = []
    for target in data:
        radec = (target[1], target[2])
        if radec not in radecs_in and radec not in outliers:
            outliers.append(radec)

    # Get XY outliers
    filename = 'output/' + cluster + '/' + date + '/phot.dat'
    data = np.loadtxt(filename)

    w = GetWCS(cluster, date)

    xmu, ymu, xystd = GetXYParams(cluster, date)

    for target in data:
        r = np.sqrt((target[0] - xmu)**2 + (target[1] - ymu)**2)
        if r > 3 * xystd:
            radec = w.all_pix2world(target[0], target[1], 0)
            radec = tuple([float(i) for i in radec])

            if radec not in outliers:
                outliers.append(radec)

    return outliers


def PMFit(cluster, date):
    data = GetGAIAInfo(cluster, date)
    outliers = GetDistanceOutliers(cluster, date)

    for target in reversed(data):
        if [target[1], target[2]] in outliers:
            data.remove(target)

    w = GetWCS(cluster, date)

    xy = [w.all_world2pix(x[1], x[2], 0) for x in data]
    x = [i[0] for i in xy]
    y = [i[1] for i in xy]
    pmra = [x[3] for x in data if not isnan(x[3])]   # mas/year
    pmdec = [x[4] for x in data if not isnan(x[4])]   # mas/year

    binning = Binning(cluster, date)

    # Plot proper motion vectors
    plt.style.use('researchpaper')
    fig, ax = plt.subplots()

    ax.quiver(x, y, pmra, pmdec, headlength=4, headaxislength=3)
    ax.set_xlim([0, 4096 / binning])
    ax.set_ylim([0, 4096 / binning])

    ax.set_xlabel('X')
    ax.set_ylabel('Y')

    ax.xaxis.set_major_locator(MultipleLocator(500 / binning))
    ax.xaxis.set_major_formatter(FormatStrFormatter('%d'))
    ax.xaxis.set_minor_locator(MultipleLocator(250 / binning))

    ax.yaxis.set_major_locator(MultipleLocator(500 / binning))
    ax.yaxis.set_major_formatter(FormatStrFormatter('%d'))
    ax.yaxis.set_minor_locator(MultipleLocator(250 / binning))

    spine_lw = 4
    [ax.spines[axis].set_linewidth(spine_lw)
     for axis in ['top', 'bottom', 'left', 'right']]

    fig.savefig('output/' + cluster + '/' + date + '/pm_vector.png')

    plt.close('all')

    # Plot RA and Dec fits
    plt.style.use('researchpaper')
    fig = plt.figure(figsize=(16, 9))

    gs = GridSpec(1, 2)

    ax1 = fig.add_subplot(gs[0, 0])
    ax2 = fig.add_subplot(gs[0, 1])

    bins = 72
    h_pmra, xpmra, _ = ax1.hist(pmra, bins=bins, color='#3f3f3f')
    h_pmdec, xpmdec, _ = ax2.hist(pmdec, bins=bins, color='#3f3f3f')
    xpmra = (xpmra[1:] + xpmra[:-1]) / 2
    xpmdec = (xpmdec[1:] + xpmdec[:-1]) / 2

    ax1.set_xlabel('PM RA')
    ax1.set_ylabel('Freq')
    ax2.set_xlabel('PM Dec')
    ax2.set_ylabel('Freq')

    median_ra = np.median(pmra)
    median_dec = np.median(pmdec)
    ra_params = hist_fit(pmra, p0=[median_ra, 1], hist_bins=72)
    dec_params = hist_fit(pmdec, p0=[median_dec, 1], hist_bins=72)

    params = [ra_params, dec_params]
    filename = 'output/' + cluster + '/' + date + '/pm_params.dat'
    with open(filename, 'w') as F:
        np.savetxt(F, params, fmt='%.3f')

    ax1.axvline(x=ra_params[0], linestyle='dashed', linewidth=1.5, color='#54ff79')
    ax1.axvline(x=ra_params[0] + 3 * ra_params[1], linestyle='dashed', linewidth=1.5,
                color='#5ec1ff')
    ax1.axvline(x=ra_params[0] - 3 * ra_params[1], linestyle='dashed', linewidth=1.5,
                color='#5ec1ff')

    ax2.axvline(x=dec_params[0], linestyle='dashed', linewidth=1.5, color='#54ff79')
    ax2.axvline(x=dec_params[0] + 3 * dec_params[1], linestyle='dashed', linewidth=1.5,
                color='#5ec1ff')
    ax2.axvline(x=dec_params[0] - 3 * dec_params[1], linestyle='dashed', linewidth=1.5,
                color='#5ec1ff')

    ax1.plot(xpmra, gauss(xpmra, *ra_params), color='#ff5454', lw=4)
    ax2.plot(xpmdec, gauss(xpmdec, *dec_params), color='#ff5454', lw=4)

    fig.savefig('output/' + cluster + '/' + date + '/pm_radec_fits.png')

    plt.close('all')


def GetPMParams(cluster, date):
    filename = 'output/' + cluster + '/' + date + '/pm_params.dat'
    try:
        params = np.loadtxt(filename)
    except IOError:
        print("'%s' does not exist." % filename)
        return

    return params.tolist()


def GetGaiaOutliers(cluster, date):
    # Get proper motion outliers
    outliers = GetDistanceOutliers(cluster, date)
    data = GetGAIAInfo(cluster, date)

    ra_params, dec_params = GetPMParams(cluster, date)
    ra_mu, ra_std, ra_a = ra_params
    dec_mu, dec_std, dec_a = dec_params

    # Get RA/Decs of stars possibly in the cluster
    radecs_in = []
    data_in = [x for x in data if abs(x[3] - ra_mu) < 3 * ra_std and
               abs(x[4] - dec_mu) < 3 * dec_std]

    for target in data_in:
        radec = (target[1], target[2])
        if radec not in radecs_in:
            radecs_in.append(radec)

    # Get RA/Decs of stars with no possibility of being in cluster
    for target in data:
        radec = (target[1], target[2])
        if radec not in radecs_in and radec not in outliers:
            outliers.append(radec)

    return outliers


def SeparatePhotometry(cluster, date):
    filename = 'output/' + cluster + '/' + date + '/phot.dat'
    data = np.loadtxt(filename).tolist()

    accepted = []
    rejected = []

    outliers = GetGaiaOutliers(cluster, date)
    w = GetWCS(cluster, date)

    for target in data:
        radec = w.all_pix2world(target[0], target[1], 0)
        radec = tuple([float(i) for i in radec])

        if radec not in outliers:
            accepted.append(target)
        else:
            rejected.append(target)

    filename = 'output/' + cluster + '/' + date + \
               '/phot_accepted.dat'
    with open(filename, 'w') as F:
        np.savetxt(F, accepted, fmt='%.3f')

    filename = 'output/' + cluster + '/' + date + \
               '/phot_rejected.dat'
    with open(filename, 'w') as F:
        np.savetxt(F, rejected, fmt='%.3f')
