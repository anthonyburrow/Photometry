from PyQt4 import QtGui, QtCore
from match import Match
from be_filter import BeFilter
from plot import Plot
from scale import Scale


class Application(QtGui.QMainWindow):

    def __init__(self):
        super(Application, self).__init__()

        self.process_type = "Single"
        self.phot_type = "psf"
        self.date = None
        self.cluster = None
        self.cooTol = 5
        self.magTol = 0.5
        self.threshold = -3.75

        # Configure window
        # self.setGeometry(50, 50, 600, 400)
        self.setWindowTitle("Photometry")
        # self.setWindowIcon(QtGui.QIcon('logo'))

        # Configure grid layouts
        centralWidget = QtGui.QWidget()
        self.mainGrid = QtGui.QGridLayout()
        centralWidget.setLayout(self.mainGrid)
        self.setCentralWidget(centralWidget)

        # Process Options
        self.photType = QtGui.QComboBox(self)
        self.photType.addItem("PSF + Aperture")
        self.photType.addItem("Aperture")
        self.photType.resize(self.photType.sizeHint())
        self.photType.activated[str].connect(self.PhotTypeChange)
        self.mainGrid.addWidget(self.photType, 0, 0)

        self.processType = QtGui.QComboBox(self)
        self.processType.addItem("Single")
        self.processType.addItem("Full")
        self.processType.resize(self.processType.sizeHint())
        self.processType.activated[str].connect(self.ProcessTypeChange)
        self.mainGrid.addWidget(self.processType, 0, 1)

        self.autoThresholdCheck = QtGui.QCheckBox("Auto Threshold", self)
        self.autoThresholdCheck.resize(self.autoThresholdCheck.sizeHint())
        self.autoThresholdCheck.toggle()
        self.autoThresholdCheck.stateChanged.connect(self.AutoThresholdCheckChange)
        self.mainGrid.addWidget(self.autoThresholdCheck, 0, 2)

        self.showCandidatesCheck = QtGui.QCheckBox("Show Candidates", self)
        self.showCandidatesCheck.resize(self.showCandidatesCheck.sizeHint())
        self.showCandidatesCheck.toggle()
        self.mainGrid.addWidget(self.showCandidatesCheck, 0, 3)

        self.lowErrorCheck = QtGui.QCheckBox("Low Error", self)
        self.lowErrorCheck.resize(self.lowErrorCheck.sizeHint())
        self.lowErrorCheck.toggle()
        self.mainGrid.addWidget(self.lowErrorCheck, 0, 4)

        # Process Functions
        self.matchCheck = QtGui.QCheckBox("Match", self)
        self.matchCheck.resize(self.matchCheck.sizeHint())
        self.mainGrid.addWidget(self.matchCheck, 1, 0)

        self.befilterCheck = QtGui.QCheckBox("Filter Be Candidates", self)
        self.befilterCheck.resize(self.befilterCheck.sizeHint())
        self.mainGrid.addWidget(self.befilterCheck, 1, 1)

        self.plotCheck = QtGui.QCheckBox("Plot", self)
        self.plotCheck.resize(self.plotCheck.sizeHint())
        self.plotCheck.stateChanged.connect(self.PlotCheckChange)
        self.mainGrid.addWidget(self.plotCheck, 1, 2)

        self.plotCMDCheck = QtGui.QCheckBox("CMD", self)
        self.plotCMDCheck.resize(self.plotCMDCheck.sizeHint())
        self.plotCMDCheck.setEnabled(False)
        self.mainGrid.addWidget(self.plotCMDCheck, 1, 3)

        self.plot2CDCheck = QtGui.QCheckBox("2CD", self)
        self.plot2CDCheck.resize(self.plot2CDCheck.sizeHint())
        self.plot2CDCheck.setEnabled(False)
        self.mainGrid.addWidget(self.plot2CDCheck, 1, 4)

        self.scaleCheck = QtGui.QCheckBox("Scale", self)
        self.scaleCheck.resize(self.scaleCheck.sizeHint())
        self.mainGrid.addWidget(self.scaleCheck, 2, 0)

        # Input specifications
        self.singleProcessDateLabel = QtGui.QLabel(self)
        self.singleProcessDateLabel.resize(self.singleProcessDateLabel.sizeHint())
        self.singleProcessDateLabel.setText("Date")
        self.mainGrid.addWidget(self.singleProcessDateLabel, 3, 0)

        self.singleProcessDate = QtGui.QLineEdit(self)
        self.singleProcessDate.setMaxLength(8)
        self.singleProcessDate.setPlaceholderText("MMDDYYYY")
        self.singleProcessDate.textChanged.connect(self.SingleProcessDateChange)
        self.mainGrid.addWidget(self.singleProcessDate, 3, 1)

        self.singleProcessClusterLabel = QtGui.QLabel(self)
        self.singleProcessClusterLabel.resize(self.singleProcessClusterLabel.sizeHint())
        self.singleProcessClusterLabel.setText("Cluster")
        self.mainGrid.addWidget(self.singleProcessClusterLabel, 4, 0)

        self.singleProcessCluster = QtGui.QLineEdit(self)
        self.singleProcessCluster.setPlaceholderText("Name")
        self.singleProcessCluster.textChanged.connect(self.SingleProcessClusterChange)
        self.mainGrid.addWidget(self.singleProcessCluster, 4, 1)

        self.cooTolInputLabel = QtGui.QLabel(self)
        self.cooTolInputLabel.resize(self.cooTolInputLabel.sizeHint())
        self.cooTolInputLabel.setText("Coo. Tol.")
        self.mainGrid.addWidget(self.cooTolInputLabel, 5, 0)

        self.cooTolInput = QtGui.QLineEdit(self)
        self.cooTolInput.setPlaceholderText("Tolerance")
        self.cooTolInput.setText(str(self.cooTol))
        self.cooTolInput.resize(self.cooTolInput.sizeHint())
        self.cooTolInput.textChanged.connect(self.CooTolChange)
        self.mainGrid.addWidget(self.cooTolInput, 5, 1)

        self.magTolInputLabel = QtGui.QLabel(self)
        self.magTolInputLabel.resize(self.magTolInputLabel.sizeHint())
        self.magTolInputLabel.setText("Mag. Tol.")
        self.mainGrid.addWidget(self.magTolInputLabel, 6, 0)

        self.magTolInput = QtGui.QLineEdit(self)
        self.magTolInput.setPlaceholderText("Tolerance")
        self.magTolInput.setText(str(self.magTol))
        self.magTolInput.resize(self.magTolInput.sizeHint())
        self.magTolInput.textChanged.connect(self.MagTolChange)
        self.mainGrid.addWidget(self.magTolInput, 6, 1)

        self.thresholdInputLabel = QtGui.QLabel(self)
        self.thresholdInputLabel.resize(self.thresholdInputLabel.sizeHint())
        self.thresholdInputLabel.setText("Threshold")
        self.mainGrid.addWidget(self.thresholdInputLabel, 7, 0)
        self.thresholdInputLabel.setEnabled(False)

        self.thresholdInput = QtGui.QLineEdit(self)
        self.thresholdInput.setPlaceholderText("Threshold")
        self.thresholdInput.setText(str(self.threshold))
        self.thresholdInput.resize(self.thresholdInput.sizeHint())
        self.thresholdInput.textChanged.connect(self.ThresholdInputChange)
        self.mainGrid.addWidget(self.thresholdInput, 7, 1)
        self.thresholdInput.setEnabled(False)

        # Process button
        self.processButton = QtGui.QPushButton("Process", self)
        self.processButton.resize(self.processButton.sizeHint())
        self.processButton.clicked.connect(self.ProcessControl)
        self.mainGrid.addWidget(self.processButton, 7, 2, 1, 3)

        self.logOutput = QtGui.QTextEdit(self)
        self.logOutput.setReadOnly(True)
        self.mainGrid.addWidget(self.logOutput, 2, 2, 5, 3)

        self.show()

    def ProcessTypeChange(self, text):
        self.process_type = text

        if text == "Single":
            self.singleProcessDate.setEnabled(True)
            self.singleProcessDate.setText(self.date)
            self.singleProcessDateLabel.setEnabled(True)
            self.singleProcessCluster.setEnabled(True)
            self.singleProcessCluster.setText(self.cluster)
            self.singleProcessClusterLabel.setEnabled(True)
        elif text == "Full":
            self.singleProcessDate.setEnabled(False)
            self.singleProcessDateLabel.setEnabled(False)
            self.singleProcessCluster.setEnabled(False)
            self.singleProcessClusterLabel.setEnabled(False)

    def PhotTypeChange(self, text):
        if text == "PSF + Aperture":
            self.phot_type = "psf"
        elif text == "Aperture":
            self.phot_type = "aperture"

    def SingleProcessDateChange(self, text):
        self.date = text

    def SingleProcessClusterChange(self, text):
        self.cluster = text

    def AutoThresholdCheckChange(self, state):
        if state == QtCore.Qt.Checked:
            self.thresholdInputLabel.setEnabled(False)
            self.thresholdInput.setEnabled(False)
        else:
            self.thresholdInputLabel.setEnabled(True)
            self.thresholdInput.setEnabled(True)

    def CooTolChange(self, text):
        try:
            self.cooTol = float(text)
        except Exception:
            pass

    def MagTolChange(self, text):
        try:
            self.magTol = float(text)
        except Exception:
            pass

    def ThresholdInputChange(self, text):
        try:
            self.threshold = float(text)
        except Exception:
            pass

    def PlotCheckChange(self, state):
        if state == QtCore.Qt.Checked:
            self.plotCMDCheck.setEnabled(True)
            self.plot2CDCheck.setEnabled(True)
        else:
            self.plotCMDCheck.setEnabled(False)
            self.plot2CDCheck.setEnabled(False)

    def ProcessControl(self):
        option = str(self.processType.currentText())

        if option == "Single":
            self.ProcessDate(self.cluster, self.date)

        elif option == "Full":
            with open("../photometry/obs_clusters.txt") as F:
                for cluster in F:
                    self.cluster = cluster
                    with open("../photometry/" + self.cluster + "/obs_dates.txt") as G:
                        for date in G:
                            self.date = date
                            self.ProcessDate(self)

                    self.ProcessCluster(self)

    def ProcessDate(self):

        input_directory = "../photometry/" + self.cluster + "/" + self.date + "/"
        output_directory = "../output/" + self.cluster + "/" + self.date + "/"

        if self.matchCheck.isChecked():
            print("Compiling all data for " + self.cluster + " on " + self.date + "...")
            match = Match(output_directory, input_directory, self.phot_type, self.cooTol, self.magTol)
            if self.lowErrorCheck.isChecked():
                match.LowError()
            else:
                match.ByExposure()

        if self.befilterCheck.isChecked():
            print("Extracting Be candidate data...")
            beFilter = BeFilter(output_directory, self.phot_type, self.autoThresholdCheck.isChecked())
            if self.lowErrorCheck.isChecked():
                beFilter.LowError()
            else:
                beFilter.Full()

        if self.plotCheck.isChecked():
            print("Generating plots...")
            plot = Plot(output_directory, self.showCandidatesCheck.isChecked(), self.lowErrorCheck.isChecked())
            if self.plotCMDCheck.isChecked():
                plot.ColorMagnitudeDiagram()
            if self.plot2CDCheck.isChecked():
                plot.TwoColorDiagram()

    def ProcessCluster(self):
        if self.scaleCheck.isChecked():
            print("Scaling observation nights for " + self.cluster)
            scale = Scale(self.cluster, self.phot_type, 10)
