import numpy as np
import matplotlib.pyplot as plt


class LowError:

    def __init__(self, cluster, date, app):
        self.cluster = cluster
        self.date = date
        self.app = app

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
        filename = "../output/" + self.cluster + "/" + self.date + "/phot_" + self.app.phot_type + ".dat"
        self.data = np.loadtxt(filename)
        V = self.data[:, 4]
        Verr = self.data[:, 5]

        # Calculate max error by std. of error
        print("\nCalculating the max error for low-error data...")

        std = np.sqrt(1 / (len(Verr) - 1)) * np.sqrt(np.sum(Verr**2))
        self.maxError = np.sqrt(2) * std   # Because max error will be for R-H and B-V

        print("  Max error found to be: ", "%.3f" % self.maxError)

        # Plot data
        self.Plot(V, Verr, std)

    def Plot(self, x, y, std):
        plt.style.use('ggplot')

        # Plot main data
        plt.plot(x, y, 'o', color='#3f3f3f', markersize=1)
        plt.title("V Err vs. V")
        plt.xlabel("V")
        plt.ylabel("V Err")
        plt.xlim([max(x), min(x)])
        plt.ylim([0, max(y) + 0.01])

        plt.hlines(std, min(x), max(x), linestyles='dashed', label='Standard Error')

        plt.legend()

        # Output
        filename = "../output/" + self.cluster + "/" + self.date + "/plots/magErr_vs_mag.png"
        plt.savefig(filename, dpi=300)

        plt.clf()
