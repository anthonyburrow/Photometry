import os.path

from .observations import ListClusters, ListDates
from .match import ProcessMatch
from .be_filter import ProcessBeFilter
from .plot import ProcessPlot
from .scale import ProcessScale, Rescale
from .gaia import ProcessGaia
from .analysis import Analysis


def Process(app):
    """Controls which dates and clusters to process.

    Depending on the process type, the application either processes
    a single date and cluster (Single), all nights for a single
    cluster, or all clusters and all dates.

    Args:
        app (Application): The GUI application object that controls processing.

    """
    option = str(app.process_type)

    if option == 'Single Date':
        SingleCluster_SingleDate(app)
    elif option == 'Single Cluster':
        SingleCluster_AllDates(app.cluster, app)
    elif option == 'Full':
        AllClusters_AllDates(app)

    print("\nComplete.")


def SingleCluster_SingleDate(app):
    """Controls processing for a single date for a single cluster.

    Args:
        app (Application): The GUI application object that controls processing.

    """
    if app.matchCheck.isChecked():
        _ProcessMatch(app.cluster, app.date, app)
    if app.befilterCheck.isChecked():
        _ProcessBeFilter(app.cluster, app.date, app, False)
    if app.plotCheck.isChecked():
        _ProcessPlot(app.cluster, app.date, app)
    if app.distanceCheck.isChecked():
        ProcessGaia(app.cluster)
    if app.summaryCheck.isChecked():
        analysis = Analysis(app.cluster, app)
        analysis.ProcessAnalysis()


def SingleCluster_AllDates(cluster, app):
    """Controls processing for every date for a single cluster.

    Args:
        app (Application): The GUI application object that controls processing.

    """
    dates = ListDates(cluster)
    baseDate = dates[0]

    # Create photometry
    if app.matchCheck.isChecked():
        for date in dates:
            print("Creating finalized photometry for %s...\n" % date)
            _ProcessMatch(cluster, date, app)

    # Scale photometry
    if app.befilterCheck.isChecked():
        for date in dates:
            print("Extracting primary Be candidates for %s on %s...\n" %
                  (cluster, date))
            _ProcessBeFilter(cluster, date, app, False)

    if app.scaleCheck.isChecked():
        for date in dates:
            print("Scaling data with primary scaling for %s on %s using " %
                  (cluster, date) + "reference %s...\n" % baseDate)
            _ProcessScale(cluster, date, app, baseDate)
        Rescale(cluster, app)

    # Analyze newly scaled photometry
    if app.distanceCheck.isChecked():
        print("Calculating cluster membership parameters...\n")
        ProcessGaia(cluster)

    if app.befilterCheck.isChecked():
        for date in dates:
            print("Extracting final Be candidates for %s on %s...\n" %
                  (cluster, date))
            _ProcessBeFilter(cluster, date, app, True)

    if app.plotCheck.isChecked():
        for date in dates:
            print("Generating plots for %s on %s..." % (cluster, date))
            _ProcessPlot(cluster, date, app)

    if app.summaryCheck.isChecked():
        print("\nCompiling Be lists and summary files...\n")
        analysis = Analysis(app.cluster, app)
        analysis.ProcessAnalysis()


def AllClusters_AllDates(app):
    """Controls processing for every date for every cluster.

    Args:
        app (Application): The GUI application object that controls processing.

    """
    clusters = ListClusters()
    for cluster in clusters:
        SingleCluster_AllDates(cluster, app)


def _ProcessMatch(cluster, date, app):
    """Processes data through matching scripts."""
    filepath = 'output/' + cluster + '/' + date
    if not os.path.exists(filepath):
        os.makedirs(filepath)

    ProcessMatch(cluster, date, app)


def _ProcessBeFilter(cluster, date, app, scaled):
    """Processes data through Be candidate filtering scripts."""
    filepath = 'output/' + cluster + '/' + date
    if not os.path.exists(filepath):
        os.makedirs(filepath)

    ProcessBeFilter(cluster, date, app, scaled)


def _ProcessPlot(cluster, date, app):
    """Processes data through plotting scripts."""
    filepath = 'output/' + cluster + '/' + date + '/plots/'
    if not os.path.exists(filepath):
        os.makedirs(filepath)

    ProcessPlot(cluster, date, app)


def _ProcessScale(cluster, date, app, baseDate):
    """Processes data through scaling scripts."""
    filepath = 'output/' + cluster + '/' + date
    if not os.path.exists(filepath):
        os.makedirs(filepath)

    ProcessScale(cluster, date, app, baseDate)
