import numpy as np
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

    def __init__(self, cluster, date, showCandidates=False, lowError=False):
        self.cluster = cluster
        self.date = date
        self.showCandidates = showCandidates
        self.lowError = lowError

        self.data = self.SetData()
        self.filtered_data = self.SetFilteredData()

        self.B = self.data[:, 2]
        self.Berr = self.data[:, 3]
        self.V = self.data[:, 4]
        self.Verr = self.data[:, 5]
        self.R = self.data[:, 6]
        self.Rerr = self.data[:, 7]
        self.H = self.data[:, 8]
        self.Herr = self.data[:, 9]

    def SetData(self):
        """Reads finalized photometry data to be plotted.

        Converts read data file to a multidimensional array used to plot data.

        Returns:
                2-dimensional array consisting of X- and Y- image coordinates, magnitudes,
                and magnitude errors.

        """
        if self.lowError:
            data = np.loadtxt("../output/" + self.cluster + "/" + self.date + "/phot_lowError")
        else:
            data = np.loadtxt("../output/" + self.cluster + "/" + self.date + "/phot.dat")

        return data

    def SetFilteredData(self):
        """Reads Be candidate photometry data to be plotted.

        Converts read data file to a multidimensional array used to plot data.

        Returns:
                2-dimensional array consisting of X- and Y- image coordinates, magnitudes,
                and magnitude errors.

        """
        if self.showCandidates:
            if self.lowError:
                filtered_data = np.loadtxt("../output/" + self.cluster + "/" + self.date + "/beList_lowError")
            else:
                filtered_data = np.loadtxt("../output/" + self.cluster + "/" + self.date + "/beList.dat")
        else:
            filtered_data = None

        return filtered_data

    def ColorMagnitudeDiagram(self):
        """Specifies the plotting of a color magnitude diagram.

        Calls to plot a color magnitude diagram with V vs. B-V axes.

        """
        print("	Creating color-magnitude diagram...")

        B_V = self.B - self.V
        B_Verr = np.sqrt(self.Berr**2 + self.Verr**2)

        self.SinglePlot(B_V, self.V, "V vs. B-V", "B-V", "V", B_Verr, self.Verr, "CMD.dat")

    def TwoColorDiagram(self):
        """Specifies the plotting of a two-color diagram.

        Calls to plot a color magnitude diagram with R-H vs. B-V axes.

        """
        print("	Creating two-color diagram...")

        B_V = self.B - self.V
        B_Verr = np.sqrt(self.Berr**2 + self.Verr**2)
        R_H = self.R - self.H
        R_Herr = np.sqrt(self.Rerr**2 + self.Herr**2)

        self.SinglePlot(B_V, R_H, "R-Ha vs. B-V", "B-V", "R-Ha", B_Verr, R_Herr, "2CD.dat")

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
        fig.savefig("../output/" + self.cluster + "/" + self.date + '/plots/' + output)
