import numpy as np
import matplotlib.pyplot as plt


class Plot:
    """Creates a Plot object, used to plot data in different ways.

    Given a data and cluster (in the form of the output directory), this class
    plots related and processed data in multiple fashions.

    Attributes:
            output_directory: Directory desired for matched output.
            showCandidates: Specify whether or not candidates are marked
            lowError: If true, only targets with lower error are shown
    """

    def __init__(self, output_directory, showCandidates=False, lowError=False):

        self.showCandidates = showCandidates
        self.lowError = lowError
        self.output_directory = output_directory

        self.data = SetData()
        self.filtered_data = SetFilteredData()

        self.B = data[:, 2]
        self.Berr = data[:, 3]
        self.V = data[:, 4]
        self.Verr = data[:, 5]
        self.R = data[:, 6]
        self.Rerr = data[:, 7]
        self.H = data[:, 8]
        self.Herr = data[:, 9]

    def SetData():
        """Reads finalized photometry data to be plotted.

        Converts read data file to a multidimensional array used to plot data.

        Returns:
                2-dimensional array consisting of X- and Y- image coordinates, magnitudes,
                and magnitude errors.

        """
        if lowError:
            data = np.loadtxt(output_directory + "phot_lowError")
        else:
            data = np.loadtxt(output_directory + "phot.dat")

        return data

    def SetFilteredData():
        """Reads Be candidate photometry data to be plotted.

        Converts read data file to a multidimensional array used to plot data.

        Returns:
                2-dimensional array consisting of X- and Y- image coordinates, magnitudes,
                and magnitude errors.

        """
        if showCandidates:
            if lowError:
                filtered_data = np.loadtxt(output_directory + "beList_lowError")
            else:
                filtered_data = np.loadtxt(output_directory + "beList.dat")
        else:
            filtered_data = None

        return filtered_data

    def ColorMagnitudeDiagram():
        """Specifies the plotting of a color magnitude diagram.

        Calls to plot a color magnitude diagram with V vs. B-V axes.

        """
        print("	Creating color-magnitude diagram...")

        B_V = B - V
        B_Verr = np.sqrt(Berr**2 + Verr**2)

        SinglePlot(B_V, V, "V vs. B-V", "B-V", "V", B_Verr, Verr, "CMD.dat")

    def TwoColorDiagram():
        """Specifies the plotting of a two-color diagram.

        Calls to plot a color magnitude diagram with R-H vs. B-V axes.

        """
        print("	Creating two-color diagram...")

        B_V = B - V
        B_Verr = np.sqrt(Berr**2 + Verr**2)
        R_H = R - H
        R_Herr = np.sqrt(Rerr**2 + Herr**2)

        SinglePlot(B_V, R_H, "R-Ha vs. B-V", "B-V", "R-Ha", B_Verr, R_Herr, "2CD.dat")

    def SinglePlot(x, y, title, x_label, y_label, x_err, y_err, output):
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
        fig.savefig(output_directory + 'plots/' + output)
