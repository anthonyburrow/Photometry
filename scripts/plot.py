import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator, FormatStrFormatter


class Plot:
    """Creates a Plot object, used to plot data in different ways.

    Given a data and cluster (in the form of the output directory), this class
    plots related and processed data in multiple fashions.

    Attributes:
            cluster:
            date:
            showCandidates: Specify whether or not candidates are marked
            lowError: If true, only targets with lower error are shown
    """

    def __init__(self, cluster, date, app):
        self.cluster = cluster
        self.date = date
        self.app = app

        self.SetData()

    def SetData(self):
        """Reads finalized photometry data to be plotted.

        Converts read data file to a multidimensional array used to plot data.
        """
        try:
            filename = "../output/" + self.cluster + "/" + self.date + "/phot_" + self.app.phot_type + "_scaled.dat"
            self.data = np.loadtxt(filename, ndmin=2)

            filename = "../output/" + self.cluster + "/" + self.date + "/phot_" + self.app.phot_type + "_lowError.dat"
            self.data_lowError = np.loadtxt(filename, ndmin=2)

            filename = "../output/" + self.cluster + "/" + self.date + "/beList_" + self.app.phot_type + ".dat"
            self.filtered_data = np.loadtxt(filename, ndmin=2)

            filename = "../output/" + self.cluster + "/" + self.date + "/beList_" + self.app.phot_type + "_lowError.dat"
            self.filtered_data_lowError = np.loadtxt(filename, ndmin=2)
        except IOError:
            print("\nFile does not exist:\n" + filename)
            return

    def ColorMagnitudeDiagram(self):
        """Specifies the plotting of a color magnitude diagram.

        Calls to plot a color magnitude diagram with V vs. B-V axes.

        """
        print("  Creating color-magnitude diagram...")

        # Plot full data
        B_V = self.data[:, 2] - self.data[:, 4]
        B_Verr = np.sqrt(self.data[:, 3]**2 + self.data[:, 5]**2)

        try:
            filtered_B_V = self.filtered_data[:, 2] - self.filtered_data[:, 4]
            filtered_V = self.filtered_data[:, 4]
        except IndexError:
            print("\nError: No data in ../output/" + self.cluster + "/" + self.date + "/beList_" + self.app.phot_type + "_lowError.dat")
            filtered_B_V = []
            filtered_V = []

        self.SinglePlot(B_V, self.data[:, 4], "V vs. B-V", "B-V", "V", B_Verr, self.data[:, 5], filtered_B_V, filtered_V, "CMD_" + self.app.phot_type)

        # Plot low-error data
        B_V_lowError = self.data_lowError[:, 2] - self.data_lowError[:, 4]
        B_Verr_lowError = np.sqrt(self.data_lowError[:, 3]**2 + self.data_lowError[:, 5]**2)

        try:
            filtered_B_V_lowError = self.filtered_data_lowError[:, 2] - self.filtered_data_lowError[:, 4]
            filtered_V_lowError = self.filtered_data_lowError[:, 4]
        except IndexError:
            print("\nError: No data in ../output/" + self.cluster + "/" + self.date + "/beList_" + self.app.phot_type + "_lowError.dat")
            filtered_B_V_lowError = []
            filtered_V_lowError = []

        self.SinglePlot(B_V_lowError, self.data_lowError[:, 4], "V vs. B-V (Low Error)", "B-V", "V", B_Verr_lowError, self.data_lowError[:, 5], filtered_B_V_lowError, filtered_V_lowError, "CMD_" + self.app.phot_type + "_lowError")

    def TwoColorDiagram(self):
        """Specifies the plotting of a two-color diagram.

        Calls to plot a color magnitude diagram with R-H vs. B-V axes.

        """
        print("  Creating two-color diagram...")

        # Plot full data
        B_V = self.data[:, 2] - self.data[:, 4]
        B_Verr = np.sqrt(self.data[:, 3]**2 + self.data[:, 5]**2)
        R_H = self.data[:, 6] - self.data[:, 8]
        R_Herr = np.sqrt(self.data[:, 7]**2 + self.data[:, 9]**2)

        try:
            filtered_B_V = self.filtered_data[:, 2] - self.filtered_data[:, 4]
            filtered_R_H = self.filtered_data[:, 6] - self.filtered_data[:, 8]
        except IndexError:
            print("\nError: No data in ../output/" + self.cluster + "/" + self.date + "/beList_" + self.app.phot_type + "_lowError.dat")
            filtered_B_V = []
            filtered_R_H = []

        self.SinglePlot(B_V, R_H, "R-Halpha vs. B-V", "B-V", "R-Halpha", B_Verr, R_Herr, filtered_B_V, filtered_R_H, "2CD_" + self.app.phot_type)

        # Plot low-error data
        B_V_lowError = self.data_lowError[:, 2] - self.data_lowError[:, 4]
        B_Verr_lowError = np.sqrt(self.data_lowError[:, 3]**2 + self.data_lowError[:, 5]**2)
        R_H_lowError = self.data_lowError[:, 6] - self.data_lowError[:, 8]
        R_Herr_lowError = np.sqrt(self.data_lowError[:, 7]**2 + self.data_lowError[:, 9]**2)

        try:
            filtered_B_V_lowError = self.filtered_data_lowError[:, 2] - self.filtered_data_lowError[:, 4]
            filtered_R_H_lowError = self.filtered_data_lowError[:, 6] - self.filtered_data_lowError[:, 8]
        except IndexError:
            print("\nError: No data in ../output/" + self.cluster + "/" + self.date + "/beList_" + self.app.phot_type + "_lowError.dat")
            filtered_B_V_lowError = []
            filtered_R_H_lowError = []

        self.SinglePlot(B_V_lowError, R_H_lowError, "R-Halpha vs. B-V (Low Error)", "B-V", "R-Halpha", B_Verr_lowError, R_Herr_lowError, filtered_B_V_lowError, filtered_R_H_lowError, "2CD_" + self.app.phot_type + "_lowError")

    def SinglePlot(self, x, y, title, x_label, y_label, x_err, y_err, be_x, be_y, output):
        """Creates a single plot of given data.

        General configuration of plotting style and other specifications, including data,
        title, labels, and error bars.  It is then output to a file in the output
        directory.

        Args:
                x: The x-axis array to be plotted.
                y: The y-axis array to be plotted.
                title: The title of the plot.
                x_label: The label of the x-axis.
                y_label: The label of the y-axis.
                x_err: The x-axis error information.
                y_err: The y-axis error information.
                be_x: The x-axis Be-filtered data.
                be_y:The y-axis Be-filtered data.
                output: The filename desired for the plot.

        """
        # plt.style.use('ggplot')

        plt.figure(figsize=(12, 9))

        # Plot main data
        plt.plot(x, y, 'o', color='#3f3f3f', markersize=12)
        # plt.title(title)
        plt.xlabel(x_label, fontsize=36)
        plt.ylabel(y_label, fontsize=36)

        plt.axes().xaxis.set_major_locator(MultipleLocator(1))
        plt.axes().xaxis.set_major_formatter(FormatStrFormatter('%d'))
        plt.axes().xaxis.set_minor_locator(MultipleLocator(0.25))

        plt.axes().yaxis.set_major_locator(MultipleLocator(2))
        plt.axes().yaxis.set_major_formatter(FormatStrFormatter('%d'))
        plt.axes().yaxis.set_minor_locator(MultipleLocator(0.5))

        plt.axes().tick_params('both', length=12, width=4, which='major', top=True, right=True, direction='in', pad=6, labelsize=30)
        plt.axes().tick_params('both', length=8, width=3, which='minor', top=True, right=True, direction='in')

        plt.axes().spines['top'].set_linewidth(4)
        plt.axes().spines['right'].set_linewidth(4)
        plt.axes().spines['bottom'].set_linewidth(4)
        plt.axes().spines['left'].set_linewidth(4)

        plt.xlim([self.app.B_VMin - 0.1, self.app.B_VMax + 2.5])
        plt.ylim([18.5 - self.app.A_v, 8.5 - self.app.A_v])
        plt.errorbar(x, y, xerr=x_err, yerr=y_err, fmt='none', ecolor='#50a0e5', elinewidth=7)

        # Overplot Be candidates
        plt.plot(be_x, be_y, 'x', color='#ff5151', markersize=15, markeredgewidth=5, label='Be Candidates')

        # Plot threshold line if 2CD
        if output == "2CD_" + self.app.phot_type or output == "2CD_" + self.app.phot_type + "_lowError":
            # plt.ylim([-5 - self.app.A_r, -1 - self.app.A_r])
            try:
                apCorr = np.loadtxt("../standards/" + self.date + "/aperture_corrections.dat")
            except IOError:
                pass

            plt.ylim([-6.5 + apCorr[2], -4 + apCorr[2]])

            plt.axes().yaxis.set_major_locator(MultipleLocator(1))
            plt.axes().yaxis.set_minor_locator(MultipleLocator(0.25))

            filename = "../output/" + self.cluster + "/" + self.date + "/thresholds_" + self.app.phot_type + ".dat"
            try:
                thresholds = np.loadtxt(filename)
                if self.app.threshold_type == "Constant":
                    file = thresholds[0]
                elif self.app.threshold_type == "Linear":
                    file = thresholds[1]
                slope = file[0]
                intercept = file[1]

                linex = np.array([self.app.B_VMin, self.app.B_VMax])
                liney = slope * linex + intercept
                plt.plot(linex, liney, '--', color='#ff5151', label='Be Threshold', linewidth=6)
            except IOError:
                print("\nNote: Thresholds have not been calculated or written to file yet and will not be displayed.")

        plt.legend(fontsize=28)
        plt.tight_layout()

        # Output
        filename = "../output/" + self.cluster + "/" + self.date + "/plots/" + output + ".png"
        plt.savefig(filename)

        plt.clf()
