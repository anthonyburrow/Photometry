from mypytools.math.gauss import gauss, bimodal, hist_fit_bimodal, gauss_2d
import numpy as np
from scipy.optimize import curve_fit
from astropy import wcs
from astropy.io import fits

import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

from .observations import ListDates


def ProcessDistances(cluster):
    print('\nFitting radial distribution for ' + cluster + '...')
    ProcessRadialDistances(cluster)

    print('Fitting x-y distributions for ' + cluster + '...')
    ProcessXYDistances(cluster)


def ProcessRadialDistances(cluster):
    all_distances = []
    sample_params = []

    # Fit radial distance distribution for each date
    print("\nCalculating " + cluster + " population distance...")
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
    # Fit x-y distributions for each date
    for date in ListDates(cluster):
        filename = 'output/' + cluster + '/' + date + '/phot_aperture.dat'
        data = np.loadtxt(filename).tolist()

        XYFit(cluster, date, data)


def GetGAIAInfo(cluster, date, dist_range=[0, 15]):
    filename = 'photometry/' + cluster + '/' + date + '/phot_dists.csv'
    try:
        data = np.genfromtxt(filename, skip_header=1, usecols=(10, 98, 99),
                             delimiter=',')   # parallax, ra, dec
    except IOError:
        print("Note: Data on distances not found for " + date)

    # Don't use negative parallax:
    data = np.array([x for x in data if x[0] > 0])
    for target in data:
        target[0] = 1 / target[0]   # convert parallax [mas] to distance [kpc]
    data = set(tuple(x) for x in data)
    data = [list(x) for x in data]   # list of unique data

    data = [x for x in data if dist_range[0] < x[0] < dist_range[1]]

    return data


def RFit(cluster, date, distances):
    # Histogram of distances
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
               color='#ff8484', label='Median = ' + '%.3f' % median_distance)

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

    # plt.cla()
    # plt.clf()

    return cluster_params


def RPopulationValues(cluster, all_distances, sample_params):
    k = len(sample_params)   # Number of observations (dates)

    sample_sizes = [x[0] for x in sample_params]
    sample_means = [x[1] for x in sample_params]
    sample_variances = [x[2] for x in sample_params]

    if k != 0:
        population_mean = np.mean(sample_means)

        z = 0
        for i in range(0, k):
            z += (sample_sizes[i] - 1) * sample_variances[i]
        variance = z / (sum(sample_sizes) - k)
        population_std = np.sqrt(variance)

        print("Cluster population found to be at " + '%.3f' % population_mean +
              " kpc +/- " + '%.3f' % population_std + " kpc")
    else:
        print("No distance data found for " + cluster +
              ".  Using d = 3.0 kpc +/- 1.0 kpc")
        population_mean = 3.0
        population_std = 1.0

    filename = 'output/' + cluster + '/Rdist_params.dat'
    with open(filename, 'w') as F:
        np.savetxt(F, [population_mean, population_std], fmt='%.3f')

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

    # plt.cla()
    # plt.clf()


def XYFit(cluster, date, data):
    filename = 'photometry/' + cluster + '/' + date + '/B1.fits'
    F = fits.getheader(filename)
    binning = F['XBINNING']

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

    fig.savefig('output/' + cluster + '/' + date + '/XYdist_distribution.png')

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
    filename = 'output/' + cluster + '/Rdist_params.dat'
    try:
        params = np.loadtxt(filename)
    except IOError:
        print("'" + filename + "'" + " does not exist.")
        return

    # Returns in form: mu, sigma
    return params.tolist()


def GetXYParams(cluster, date):
    filename = 'output/' + cluster + '/' + date + '/XYdist_params.dat'
    try:
        params = np.loadtxt(filename)
    except IOError:
        print("'" + filename + "'" + " does not exist.")
        return

    # Returns in form: mu, sigma
    return params.tolist()


def GetDistanceOutliers(cluster, date):
    # Get radial outliers
    R_outliers = []

    data = GetGAIAInfo(cluster, date)

    radec = [[x[1], x[2]] for x in data]
    radec = set(tuple(x) for x in radec)
    radec = [list(x) for x in radec]   # list of unique ra/dec for iteration

    rmu, rstd = GetRParams(cluster)

    for target in radec:
        d = []
        for line in [x for x in data
                     if x[1] == target[0] and x[2] == target[1]]:
            d.append(line[0])
        if not any(abs(x - rmu) < 3 * rstd for x in d):
            R_outliers.append(target)

    # Get XY outliers
    XY_outliers = []

    filename = 'output/' + cluster + '/' + date + '/phot_aperture.dat'
    data = np.loadtxt(filename)

    filename = 'photometry/' + cluster + '/' + date + '/B1_wcs.fits'
    w = wcs.WCS(filename)

    xmu, ymu, xystd = GetXYParams(cluster, date)

    for target in data:
        if np.sqrt((target[0] - xmu)**2 + (target[1] - ymu)**2) > 3 * xystd:
            radec = w.all_pix2world(target[0], target[1], 0)
            ra = float(radec[0])
            dec = float(radec[1])
            if [ra, dec] not in R_outliers:
                XY_outliers.append([ra, dec])

    return R_outliers + XY_outliers
