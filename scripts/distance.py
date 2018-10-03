import numpy as np
import matplotlib.pyplot as plt
# from matplotlib.ticker import MultipleLocator, FormatStrFormatter
from scipy.optimize import curve_fit
# import pandas as pd
from observations import ListDates


class Distance:

    def __init__(self, cluster):
        self.cluster = cluster

        self.sample_params = []
        self.all_distances = []

    def process(self):
        print("Calculating " + self.cluster + " population distance...")
        for date in ListDates(self.cluster):
            self.process_date(date)

        self.population_values()

    def process_date(self, date):
        # Get info from file given by gaia archive
        filename = '../photometry/' + self.cluster + '/' + date + '/phot_dists.csv'
        try:
            data = np.genfromtxt(filename, skip_header=1, usecols=(10, 98, 99), delimiter=',')   # parallax, ra, dec
        except IOError:
            print("  Data on distances not found for " + date)
            return

        data = np.array([x for x in data if x[0] > 0])   # don't use negative parallax
        for target in data:
            target[0] = 1 / target[0]   # convert parallax [mas] to distance [kpc]
        distances = list(set([x[0] for x in data]))   # get rid of duplicates
        dMax = 15
        distances = [x for x in distances if x < dMax]
        self.all_distances.extend(distances)

        self.fit(distances, date)

    def fit(self, distances, date):
        # Histogram of distances
        plt.figure(figsize=(12, 9))

        low_lim = 0
        up_lim = 10

        y, x, _ = plt.hist(distances, bins=75, range=[low_lim, up_lim], color="#3f3f3f")
        plt.xlabel("Distance (kpc)", fontsize=36)
        plt.ylabel("Frequency", fontsize=36)

        plt.xlim(low_lim, up_lim)

        median_distance = np.median(distances)
        plt.axvline(x=median_distance, linestyle='dashed', linewidth=1.5, color="#ff8484", label='Median = ' + '%.3f' % median_distance)

        # plt.axes().xaxis.set_major_locator(MultipleLocator(1))
        # plt.axes().xaxis.set_major_formatter(FormatStrFormatter('%.1f'))
        # plt.axes().xaxis.set_minor_locator(MultipleLocator(0.25))

        # Fit bimodal Gaussians to data
        x = (x[1:] + x[:-1]) / 2  # for len(x)==len(y)

        def gauss(x, mu, sigma, A):
            return A * np.exp(-(x - mu)**2 / 2 / sigma**2)

        def bimodal(x, mu1, sigma1, A1, mu2, sigma2, A2):
            return gauss(x, mu1, sigma1, A1) + gauss(x, mu2, sigma2, A2)

        ymax = max(y)
        # p0 = (0, 1, ymax / 2, np.median(distances), 1, ymax)
        p0 = (0.8, 1, 38, 3.0, 1, ymax)
        params, cov = curve_fit(bimodal, x, y, p0, bounds=((low_lim, 0, 0, low_lim, 0, 0), (up_lim, up_lim - low_lim, np.inf, up_lim, up_lim - low_lim, np.inf)))
        # sigma = np.sqrt(np.diag(cov))

        # Use curve closer to median (would use amplitude comparison instead of mean but it was acting wonky)
        if abs(params[0] - median_distance) < abs(params[3] - median_distance):
            cluster_params = params[0:3]
            curve = 'G1'
        else:
            cluster_params = params[3:6]
            curve = 'G2'

        plt.axvline(x=cluster_params[0], linestyle='dashed', linewidth=1.5, color="#54ff79", label=r'$\mu$' + '(' + curve + ') = ' + '%.3f' % cluster_params[0])
        plt.axvline(x=cluster_params[0] + 3 * cluster_params[1], linestyle='dashed', linewidth=1.5, color="#5ec1ff")
        plt.axvline(x=cluster_params[0] - 3 * cluster_params[1], linestyle='dashed', linewidth=1.5, color="#5ec1ff")

        plt.plot(x, bimodal(x, *params), color='#ff5454', lw=4)
        plt.plot(x, gauss(x, *params[0:3]), color='#5495ff', ls='--', lw=4, label='G1')
        plt.plot(x, gauss(x, *params[3:6]), color='#54ff79', ls='--', lw=4, label='G2')

        # print(pd.DataFrame(data={'params': params, 'sigma': sigma}, index=bimodal.__code__.co_varnames[1:]))

        # Record sample values for population calculation
        self.sample_params.append([len(distances), cluster_params[0], cluster_params[1]**2, cluster_params[2]])

        # Axes specs
        plt.axes().tick_params('both', length=12, width=4, which='major', top=True, right=True, direction='in', pad=6, labelsize=30)
        plt.axes().tick_params('both', length=8, width=3, which='minor', top=True, right=True, direction='in')

        plt.axes().spines['top'].set_linewidth(4)
        plt.axes().spines['right'].set_linewidth(4)
        plt.axes().spines['bottom'].set_linewidth(4)
        plt.axes().spines['left'].set_linewidth(4)

        plt.legend(fontsize=28)
        plt.tight_layout()

        filename = '../output/' + self.cluster + '/' + date + '/distance_distribution.png'
        plt.savefig(filename)
        plt.clf()

    def population_values(self):
        k = len(self.sample_params)   # Number of observations (dates)

        sample_sizes = [x[0] for x in self.sample_params]
        sample_means = [x[1] for x in self.sample_params]
        sample_variances = [x[2] for x in self.sample_params]

        if k != 0:
            population_mean = np.mean(sample_means)

            z = 0
            for i in range(0, k):
                z += (sample_sizes[i] - 1) * sample_variances[i]
            variance = z / (sum(sample_sizes) - k)
            population_std = np.sqrt(variance)

            print("Cluster population found to be at " + '%.3f' % population_mean + " kpc +/- " + '%.3f' % population_std + " kpc")
        else:
            print("No distance data found for " + self.cluster + ".  Using d = 3.0 kpc +/- 1.0 kpc")
            population_mean = 3.0
            population_std = 1.0

        filename = '../output/' + self.cluster + '/distance_params.dat'
        with open(filename, 'w') as F:
            F.write(str(population_mean) + '\n' + str(population_std))

        # Plot population data
        plt.figure(figsize=(12, 9))

        low_lim = 0
        up_lim = 10
        bins = 100

        y, x, _ = plt.hist(self.all_distances, bins=100, range=[low_lim, up_lim], color="#3f3f3f")
        x = (x[1:] + x[:-1]) / 2  # for len(x)==len(y)
        plt.clf()

        plt.bar(x, y / k, (up_lim - low_lim) / bins, align='center', color="#3f3f3f")

        plt.xlabel("Distance (kpc)", fontsize=36)
        plt.ylabel("Frequency", fontsize=36)

        plt.xlim(low_lim, up_lim)

        def gauss(x, mu, sigma, A):
            return A * np.exp(-(x - mu)**2 / 2 / sigma**2)

        def superpose(x, params):
            g = np.zeros_like(x)
            for param in params:
                g += gauss(x, *param[1:])
            g = g / len(params)   # Normalize
            return g

        # Calculated population mean
        plt.axvline(x=population_mean, linestyle='dashed', linewidth=1.5, color="#54ff79", label=r'$\mu$' + ' = ' + '%.3f' % population_mean)
        # Calculated population std (2-sigma)
        plt.axvline(x=population_mean + 2 * population_std, linestyle='dashed', linewidth=1.5, color="#5ec1ff", label='2-sigma')
        plt.axvline(x=population_mean - 2 * population_std, linestyle='dashed', linewidth=1.5, color="#5ec1ff")
        # Calculated population std (3-sigma)
        plt.axvline(x=population_mean + 3 * population_std, linestyle='dashed', linewidth=1.5, color="#dc5eff", label='3-sigma')
        plt.axvline(x=population_mean - 3 * population_std, linestyle='dashed', linewidth=1.5, color="#dc5eff")

        # Weighted sum of gaussians
        plt.plot(x, superpose(x, self.sample_params), color='#f44242', ls='-', lw=4)

        # Axes specs
        plt.axes().tick_params('both', length=12, width=4, which='major', top=True, right=True, direction='in', pad=6, labelsize=30)
        plt.axes().tick_params('both', length=8, width=3, which='minor', top=True, right=True, direction='in')

        plt.axes().spines['top'].set_linewidth(4)
        plt.axes().spines['right'].set_linewidth(4)
        plt.axes().spines['bottom'].set_linewidth(4)
        plt.axes().spines['left'].set_linewidth(4)

        plt.legend(fontsize=28)
        plt.tight_layout()

        filename = '../output/' + self.cluster + '/distance_distribution.png'
        plt.savefig(filename)
        plt.clf()
