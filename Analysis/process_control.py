import os.path

from .observations import ListClusters, ListDates
from .match import ProcessMatch
from .low_error import ProcessLowError
from .be_filter import ProcessBeFilter
from .plot import ProcessPlot
from .scale import ProcessScale, Rescale
from .distance import ProcessDistances
from .analysis import ProcessAnalysis


def Process(app):
    """Controls which dates and clusters to process.

    Depending on the process type, the application
    either processes a single date/cluster (Single)
    or processes all nights and clusters with finalized photometry.

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
    if app.matchCheck.isChecked():
        _ProcessMatch(app.cluster, app.date, app)
    if app.lowErrorCheck.isChecked():
        _ProcessLowError(app.cluster, app.date, app, 'phot')
    if app.befilterCheck.isChecked():
        _ProcessBeFilter(app.cluster, app.date, app, False)
    if app.plotCheck.isChecked():
        _ProcessPlot(app.cluster, app.date, app)
    if app.distanceCheck.isChecked():
        ProcessDistances(app.cluster)
    if app.summaryCheck.isChecked():
        ProcessAnalysis(app.cluster, app)


def SingleCluster_AllDates(cluster, app):
    dates = ListDates(cluster)
    baseDate = dates[0]

    # Create photometry
    if app.matchCheck.isChecked():
        for date in dates:
            _ProcessMatch(cluster, date, app)

    # Scale photometry
    if app.lowErrorCheck.isChecked():
        for date in dates:
            _ProcessLowError(cluster, date, app, 'phot')

    if app.befilterCheck.isChecked():
        for date in dates:
            _ProcessBeFilter(cluster, date, app, False)

    if app.scaleCheck.isChecked():
        for date in dates:
            _ProcessScale(cluster, date, app, baseDate)
        Rescale(cluster, app)

    # Analyze newly scaled photometry
    if app.lowErrorCheck.isChecked():
        for date in dates:
            _ProcessLowError(cluster, date, app, 'phot_scaled')

    if app.befilterCheck.isChecked():
        for date in dates:
            _ProcessBeFilter(cluster, date, app, True)
    if app.lowErrorCheck.isChecked():
        for date in dates:
            _ProcessLowError(cluster, date, app, 'beList_scaled')

    if app.plotCheck.isChecked():
        for date in dates:
            filepath = 'output/' + cluster + '/' + date + '/plots/'
            if not os.path.exists(filepath):
                os.makedirs(filepath)
            _ProcessPlot(cluster, date, app)

    if app.distanceCheck.isChecked():
        ProcessDistances(cluster)

    if app.summaryCheck.isChecked():
        ProcessAnalysis(cluster, app)


def AllClusters_AllDates(app):
    """Calls each process type for each date and for each cluster."""
    clusters = ListClusters()
    for cluster in clusters:
        SingleCluster_AllDates(cluster, app)


def _ProcessMatch(cluster, date, app):
    """Processes data through the matching scripts."""
    filepath = 'output/' + cluster + '/' + date
    if not os.path.exists(filepath):
        os.makedirs(filepath)

    ProcessMatch(cluster, date, app)


def _ProcessLowError(cluster, date, app, file):
    path = 'output/' + cluster + '/' + date + '/'
    ProcessLowError(path, file, app)


def _ProcessBeFilter(cluster, date, app, scaled):
    """Processes data through the Be candidate filtering scripts."""
    filepath = 'output/' + cluster + '/' + date
    if not os.path.exists(filepath):
        os.makedirs(filepath)

    ProcessBeFilter(cluster, date, app, scaled)


def _ProcessPlot(cluster, date, app):
    """Processes data through the plotting scripts."""
    filepath = 'output/' + cluster + '/' + date + '/plots/'
    if not os.path.exists(filepath):
        os.makedirs(filepath)

    print("\nGenerating plots for " + cluster + " on " + date + "...\n")
    ProcessPlot(cluster, date, app)


def _ProcessScale(cluster, date, app, baseDate):
    """Processes data through the scaling scripts."""
    filepath = 'output/' + cluster + '/' + date
    if not os.path.exists(filepath):
        os.makedirs(filepath)

    ProcessScale(cluster, date, app, baseDate)
