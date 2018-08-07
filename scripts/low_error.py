import numpy as np
import os.path
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator, FormatStrFormatter


class LowError:

    def __init__(self, cluster, date, app, scaled):
        self.cluster = cluster
        self.date = date
        self.app = app
        self.scaled = scaled

        self.MaxError()

    def Process(self):
        lowError_data = []

        Rerr = np.array(self.data)[:, 7]
        Herr = np.array(self.data)[:, 9]

        R_Herr = np.sqrt(Rerr**2 + Herr**2)

        for i in range(0, len(self.data)):
            if R_Herr[i] < self.maxError:
                lowError_data.append(self.data[i])

        # Output to file
        filename = "../output/" + self.cluster + "/" + self.date + "/phot_" + self.app.phot_type + "_lowError.dat"
        with open(filename, 'w') as F:
            np.savetxt(F, lowError_data, fmt='%.3f')

        return lowError_data

    def MaxError(self):
        # Read data
        if self.scaled:
            filename = "../output/" + self.cluster + "/" + self.date + "/phot_" + self.app.phot_type + "_scaled.dat"
        else:
            filename = "../output/" + self.cluster + "/" + self.date + "/phot_" + self.app.phot_type + ".dat"

        self.data = np.loadtxt(filename)
        R_H = self.data[:, 6] - self.data[:, 8]
        R_Herr = np.sqrt(self.data[:, 7]**2 + self.data[:, 9]**2)

        # Calculate max error by std. of error
        print("\nCalculating the max error for low-error data...")

        std = np.sqrt(1 / (len(R_Herr) - 1)) * np.sqrt(np.sum(R_Herr**2))
        self.maxError = np.sqrt(2) * std

        print("  Max error found to be: ", "%.3f" % self.maxError)

        # Plot data
        self.Plot(R_H, R_Herr, std)

    def Plot(self, x, y, std):
        # plt.style.use('ggplot')

        plt.figure(figsize=(12, 9))

        # Plot main data
        plt.plot(x, y, 'o', color='#3f3f3f', markersize=12)
        # plt.title("R-H Err vs. R-H")
        plt.xlabel("R-H", fontsize=36)
        plt.ylabel("R-H Err", fontsize=36)
        plt.xlim([max(x), min(x)])
        plt.ylim([min(y) - 0.01, max(y) + 0.01])

        plt.axes().xaxis.set_major_locator(MultipleLocator(2))
        plt.axes().xaxis.set_major_formatter(FormatStrFormatter('%d'))
        plt.axes().xaxis.set_minor_locator(MultipleLocator(0.5))

        plt.axes().yaxis.set_major_locator(MultipleLocator(0.02))
        plt.axes().yaxis.set_major_formatter(FormatStrFormatter('%.2f'))
        plt.axes().yaxis.set_minor_locator(MultipleLocator(0.005))

        plt.axes().tick_params('both', length=12, width=4, which='major', top=True, right=True, direction='in', pad=6, labelsize=30)
        plt.axes().tick_params('both', length=8, width=3, which='minor', top=True, right=True, direction='in')

        plt.axes().spines['top'].set_linewidth(4)
        plt.axes().spines['right'].set_linewidth(4)
        plt.axes().spines['bottom'].set_linewidth(4)
        plt.axes().spines['left'].set_linewidth(4)

        plt.hlines(std, min(x), max(x), linestyles='dashed', label='Standard Error')

        plt.legend(fontsize=28)
        plt.tight_layout()

        # Output
        if not os.path.exists("../output/" + self.cluster + "/" + self.date + "/plots/"):
            os.makedirs("../output/" + self.cluster + "/" + self.date + "/plots/")

        filename = "../output/" + self.cluster + "/" + self.date + "/plots/magErr_vs_mag_" + self.app.phot_type + ".png"
        plt.savefig(filename)

        plt.clf()
