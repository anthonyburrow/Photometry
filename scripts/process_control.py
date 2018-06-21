from match import Match
from be_filter import BeFilter
from plot import Plot
from scale import Scale


class ProcessControl:

    def __init__(self, app):
        self.app = app

    # Cluster Specific
    def AllClusters(self, func, *args):
        with open("../photometry/obs_clusters.txt") as F:
            clusters = F.readlines()
            for cluster in clusters:
                self.SingleCluster(cluster, func, *args)

    def SingleCluster(cluster, func, *args):
        func(cluster, *args)

    # Date specific
    def AllClusters_AllDates(self, func, *args):
        with open("../photometry/obs_clusters.txt") as F:
            clusters = F.readlines()
            for cluster in clusters:
                self.SingleCluster_AllDates(cluster, func, *args)

    def SingleCluster_AllDates(self, cluster, func, *args):
        with open("../photometry/" + cluster + "/obs_dates.txt") as F:
            dates = F.readlines()
            for date in dates:
                self.SingleCluster_SingleDate(cluster, date, func, *args)

    def SingleCluster_SingleDate(cluster, date, func, *args):
        func(cluster, date, *args)

    def Process(self):
        """Controls which dates and clusters to process.

        Depending on the process type, the application either processes a single date/cluster (Single)
        or processes all nights and clusters with finalized photometry.

        """
        option = str(self.app.processType.currentText())

        if option == "Single":
            self.ProcessMatch(self.app.cluster, self.app.date)
            self.ProcessBeFilter(self.app.cluster, self.app.date)
            self.ProcessPlot(self.app.cluster, self.app.date)
            self.ProcessScale(self.app.cluster, self.app.date)

        elif option == "Full":
            self.AllClusters_AllDates(self.ProcessMatch)
            self.AllClusters_AllDates(self.ProcessBeFilter)
            self.AllClusters_AllDates(self.ProcessPlot)
            self.AllClusters_AllDates(self.ProcessScale)

    def ProcessMatch(self, cluster, date):
        if self.app.matchCheck.isChecked():
            print("Compiling all data for " + cluster + " on " + date + "...")
            match = Match(cluster, date, self.app.phot_type, self.app.cooTol, self.app.magTol)
            if self.app.lowErrorCheck.isChecked():
                match.LowError()
            else:
                match.ByExposure()

    def ProcessBeFilter(self, cluster, date, scaled=False):
        if self.app.befilterCheck.isChecked():
            print("Extracting Be candidates for " + cluster + " on " + date + "...")
            beFilter = BeFilter(cluster, date, self.app.phot_type, self.app.autoThresholdCheck.isChecked(), self.app.threshold_type, self.app.threshold, [self.app.B_VMin, self.app.B_VMax], scaled)
            if self.app.lowErrorCheck.isChecked():
                beFilter.LowError()
            else:
                beFilter.Full()

    def ProcessPlot(self, cluster, date):
        if self.app.plotCheck.isChecked():
            print("Generating plots for " + cluster + " on " + date + "...")
            plot = Plot(cluster, date, self.app.showCandidatesCheck.isChecked(), self.app.lowErrorCheck.isChecked())
            if self.app.plotCMDCheck.isChecked():
                plot.ColorMagnitudeDiagram()
            if self.app.plot2CDCheck.isChecked():
                plot.TwoColorDiagram()

    def ProcessScale(self, cluster, date):
        """Processes all data from a single cluster.

        For the entire cluster, each night is scaled to each other using the Scale module.

        """
        if self.app.scaleCheck.isChecked():
            # Set up scale object and base date data
            scale = Scale(cluster, date, self.app.phot_type, 10)

            # Scale only if not the reference date
            with open("../photometry/" + cluster + "/obs_dates.txt") as F:  # Establish first date as scaling base
                baseDate = F.readline()
                if date != baseDate:
                    print("Scaling data for " + cluster + " on " + date)
                    scale.Scale()

        # Filter scaled data
        self.ProcessBeFilter(cluster, date, True)
