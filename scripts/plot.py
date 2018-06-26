import numpy as np
import os.path
import matplotlib.pyplot as plt


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
        self.data = np.loadtxt("../output/" + self.cluster + "/" + self.date + "/phot_" + self.app.phot_type + ".dat")
        # self.data_lowError = np.loadtxt("../output/" + self.cluster + "/" + self.date + "/phot_" + self.app.phot_type + "_lowError.dat")
        # if self.app.showCandidatesCheck.isChecked():
        #     self.filtered_data = np.loadtxt("../output/" + self.cluster + "/" + self.date + "/beList.dat")
        #     self.filtered_data_lowError = np.loadtxt("../output/" + self.cluster + "/" + self.date + "/beList_lowError")
        # else:
        #     self.filtered_data = []
        #     self.filtered_data_lowError = []

    def ColorMagnitudeDiagram(self):
        """Specifies the plotting of a color magnitude diagram.

        Calls to plot a color magnitude diagram with V vs. B-V axes.

        """
        print(" Creating color-magnitude diagram...")

        B_V = self.data[:, 2] - self.data[:, 4]
        B_Verr = np.sqrt(self.data[:, 3]**2 + self.data[:, 5]**2)

        B_V_lowError = self.data_lowError[:, 2] - self.data_lowError[:, 4]
        B_Verr_lowError = np.sqrt(self.data_lowError[:, 3]**2 + self.data_lowError[:, 5]**2)

        self.SinglePlot(B_V, self.data[:, 4], "V vs. B-V", "B-V", "V", B_Verr, self.data[:, 5], "CMD.dat")
        self.SinglePlot(B_V_lowError, self.data_lowError[:, 4], "V vs. B-V (Low Error)", "B-V", "V", B_Verr_lowError, self.data_lowError[:, 5], "CMD_lowError.dat")

    def TwoColorDiagram(self):
        """Specifies the plotting of a two-color diagram.

        Calls to plot a color magnitude diagram with R-H vs. B-V axes.

        """
        print(" Creating two-color diagram...")

        B_V = self.data[:, 2] - self.data[:, 4]
        # B_Verr = np.sqrt(self.data[:, 3]**2 + self.data[:, 5]**2)
        R_H = self.data[:, 6] - self.data[:, 8]
        # R_Herr = np.sqrt(self.data[:, 7]**2 + self.data[:, 9]**2)

        # Plot data
        plt.style.use('ggplot')

        fig = plt.figure(1)
        ax = fig.add_subplot(111)

        ax.plot(B_V, R_H, 'o')
        ax.set_title("R-Ha vs. B-V")
        ax.set_xlim([self.app.B_VMin - 0.1, 3])
        ax.set_ylim([-5, -1])

        filename = "../output/" + self.cluster + "/" + self.date + "/thresholds.dat"
        if os.path.isfile(filename):
            thresholds = np.loadtxt(filename)
            if self.app.threshold_type == "Constant":
                file = thresholds[0]
            elif self.app.threshold_type == "Linear":
                file = thresholds[1]
            slope = file[0]
            intercept = file[1]

            x = np.array([self.app.B_VMin, self.app.B_VMax])
            y = slope * x + intercept
            plt.plot(x, y, '--')

        # plt.show()
        filename = "../output/" + self.cluster + "/" + self.date + '/plots/2CD.png'
        fig.savefig(filename, dpi=300)

        # B_V_lowError = self.data_lowError[:, 2] - self.data_lowError[:, 4]
        # B_Verr_lowError = np.sqrt(self.data_lowError[:, 3]**2 + self.data_lowError[:, 5]**2)
        # R_H_lowError = self.data_lowError[:, 6] - self.data_lowError[:, 8]
        # R_Herr_lowError = np.sqrt(self.data_lowError[:, 7]**2 + self.data_lowError[:, 9]**2)

    def SinglePlot(self, x, y, title, x_label, y_label, x_err, y_err, output):
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
                output: The filename desired for the plot.

        """
        # May need to change a lot with styling,
        # including adding lines to signify candidate limits,
        # toggling and implementing candidate visibility,
        # implementing error bars,
        # and other general changes
        plt.style.use('ggplot')

        fig = plt.figure(1)
        ax = fig.add_subplot(111)

        ax.plot(x, y, 'o')
        ax.set_title()

        # plt.show()
        filename = "../output/" + self.cluster + "/" + self.date + '/plots/' + output
        fig.savefig(filename)
