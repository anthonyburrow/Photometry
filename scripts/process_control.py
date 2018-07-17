import os.path
from observations import Observations
from low_error import LowError
from analysis import Analysis


class ProcessControl:
    """Controls which processes are to occur for a run.

    Attributes:
            app: GUI application that specifies parameters.
    """

    def __init__(self, app):
        self.app = app

    def Process(self):
        """Controls which dates and clusters to process.

        Depending on the process type, the application either processes a single date/cluster (Single)
        or processes all nights and clusters with finalized photometry.

        """
        option = str(self.app.process_type)

        if option == "Single":
            self.ProcessMatch(self.app.cluster, self.app.date)
            self.ProcessBeFilter(self.app.cluster, self.app.date)
            self.ProcessPlot(self.app.cluster, self.app.date)
            self.ProcessScale(self.app.cluster, self.app.date)
            if self.app.summaryCheck.isChecked():
                Analysis(self.app.cluster, self.app).Summary()

        elif option == "Full":
            self.AllClusters_AllDates()

        print("\nComplete.")

    def AllClusters_AllDates(self):
        """Calls each process type for each date and for each cluster."""
        clusters = Observations().ListClusters()
        for cluster in clusters:
            dates = Observations().ListDates(cluster)
            for date in dates:
                self.ProcessMatch(cluster, date)
                self.ProcessBeFilter(cluster, date)
                self.ProcessPlot(cluster, date)
                self.ProcessScale(cluster, date)
            if self.app.summaryCheck.isChecked():
                Analysis(cluster, self.app).Summary()

    def ProcessMatch(self, cluster, date):
        """Processes data through the matching scripts."""
        if not os.path.exists("../output/" + cluster + "/" + date):
            os.makedirs("../output/" + cluster + "/" + date)

        if self.app.matchCheck.isChecked():
            from match import Match

            print("\nCompiling all data for " + cluster + " on " + date + "...\n")
            match = Match(cluster, date, self.app)
            match.ByExposure()
            lowError = LowError(cluster, date, self.app)
            lowError.Process()

    def ProcessBeFilter(self, cluster, date):
        """Processes data through the Be candidate filtering scripts."""
        if not os.path.exists("../output/" + cluster + "/" + date):
            os.makedirs("../output/" + cluster + "/" + date)

        if self.app.befilterCheck.isChecked():
            from be_filter import BeFilter

            print("\nExtracting Be candidates for " + cluster + " on " + date + "...\n")
            beFilter = BeFilter(cluster, date, self.app)
            beFilter.Process()

    def ProcessPlot(self, cluster, date):
        """Processes data through the plotting scripts."""
        if not os.path.exists("../output/" + cluster + "/" + date + "/plots/"):
            os.makedirs("../output/" + cluster + "/" + date + "/plots/")

        if self.app.plotCheck.isChecked():
            from plot import Plot

            print("\nGenerating plots for " + cluster + " on " + date + "...\n")
            plot = Plot(cluster, date, self.app)
            if self.app.plotCMDCheck.isChecked():
                plot.ColorMagnitudeDiagram()
            if self.app.plot2CDCheck.isChecked():
                plot.TwoColorDiagram()

    def ProcessScale(self, cluster, date):
        """Processes data through the scaling scripts."""
        if not os.path.exists("../output/" + cluster + "/" + date):
            os.makedirs("../output/" + cluster + "/" + date)

        if self.app.scaleCheck.isChecked():
            from scale import Scale

            # Set up scale object and base date data
            scale = Scale(cluster, date, self.app)   # This creates reference night 'scaled' file

            # Scale only if not the reference date
            baseDate = Observations().ListDates(cluster)[0]   # Establish first date as scaling base
            if date != baseDate:   # Don't need to scale the reference date
                print("\nScaling data for " + cluster + " on " + date + "...\n")
                scale.Scale()
