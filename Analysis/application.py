from PyQt4 import QtGui, QtCore

from .config.config import GetSettings

from .process_control import Process


class Application(QtGui.QMainWindow):

    def __init__(self):
        super(Application, self).__init__()

        settings = GetSettings()

        self.process_type = settings['process_options']['process_type']
        self.threshold_type = settings['process_options']['threshold_type']
        self.date = settings['process_options']['date']
        self.cluster = settings['process_options']['cluster']
        self.cooTol = settings['data']['cooTol']
        self.magTol = settings['data']['magTol']
        self.threshold = settings['data']['threshold']
        self.B_VMax = settings['data']['B_VMax']
        self.B_VMin = settings['data']['B_VMin']

        # Configure window
        # self.setGeometry(50, 50, 600, 400)
        self.setWindowTitle("Photometry")
        # self.setWindowIcon(QtGui.QIcon('logo'))

        # Configure grid layouts
        centralWidget = QtGui.QWidget()
        self.mainGrid = QtGui.QGridLayout()
        centralWidget.setLayout(self.mainGrid)
        self.setCentralWidget(centralWidget)

        # Process type
        self.processType = QtGui.QComboBox(self)
        self.processType.addItem("Single Cluster")
        self.processType.addItem("Full")
        self.processType.addItem("Single Date")
        self.processType.resize(self.processType.sizeHint())
        self.processType.activated[str].connect(self.ProcessTypeChange)
        self.mainGrid.addWidget(self.processType, 0, 0)

        # Threshold type
        self.thresholdType = QtGui.QComboBox(self)
        self.thresholdType.addItem("Linear")
        self.thresholdType.addItem("Constant")
        self.thresholdType.resize(self.thresholdType.sizeHint())
        self.thresholdType.activated[str].connect(self.ThresholdTypeChange)
        self.mainGrid.addWidget(self.thresholdType, 0, 1)

        # Auto threshold checkbox
        self.autoThresholdCheck = QtGui.QCheckBox("Auto Threshold", self)
        self.autoThresholdCheck.resize(self.autoThresholdCheck.sizeHint())
        self.autoThresholdCheck.toggle()
        self.autoThresholdCheck.stateChanged.connect(
            self.AutoThresholdCheckChange)
        self.mainGrid.addWidget(self.autoThresholdCheck, 0, 2)

        # Run match process
        self.matchCheck = QtGui.QCheckBox("Match", self)
        self.matchCheck.resize(self.matchCheck.sizeHint())
        if settings['processes']['match']:
            self.matchCheck.toggle()
        self.mainGrid.addWidget(self.matchCheck, 1, 0)

        # Run Be filter process
        self.befilterCheck = QtGui.QCheckBox("Filter Be Candidates", self)
        self.befilterCheck.resize(self.befilterCheck.sizeHint())
        if settings['processes']['be_filter']:
            self.befilterCheck.toggle()
        self.mainGrid.addWidget(self.befilterCheck, 1, 1)

        # Run scale process
        self.scaleCheck = QtGui.QCheckBox("Scale", self)
        self.scaleCheck.resize(self.scaleCheck.sizeHint())
        if settings['processes']['scale']:
            self.scaleCheck.toggle()
        self.mainGrid.addWidget(self.scaleCheck, 1, 2)

        # Run plot process
        self.plotCheck = QtGui.QCheckBox("Plot", self)
        self.plotCheck.resize(self.plotCheck.sizeHint())
        if settings['processes']['plot']:
            self.plotCheck.toggle()
        self.mainGrid.addWidget(self.plotCheck, 1, 4)

        # Run distance process
        self.distanceCheck = QtGui.QCheckBox("Distances", self)
        self.distanceCheck.resize(self.distanceCheck.sizeHint())
        if settings['processes']['distance']:
            self.distanceCheck.toggle()
        self.mainGrid.addWidget(self.distanceCheck, 1, 3)

        # Run analysis process
        self.summaryCheck = QtGui.QCheckBox("Summary", self)
        self.summaryCheck.resize(self.summaryCheck.sizeHint())
        if settings['processes']['analysis']:
            self.summaryCheck.toggle()
        self.mainGrid.addWidget(self.summaryCheck, 2, 0)

        # Label for manual single date
        self.singleProcessDateLabel = QtGui.QLabel(self)
        self.singleProcessDateLabel.resize(
            self.singleProcessDateLabel.sizeHint())
        self.singleProcessDateLabel.setText("Date")
        self.mainGrid.addWidget(self.singleProcessDateLabel, 3, 0)
        self.singleProcessDateLabel.setEnabled(False)

        # Input for manual single date
        self.singleProcessDate = QtGui.QLineEdit(self)
        self.singleProcessDate.setMaxLength(8)
        self.singleProcessDate.setPlaceholderText("MMDDYYYY")
        self.singleProcessDate.textChanged.connect(
            self.SingleProcessDateChange)
        self.mainGrid.addWidget(self.singleProcessDate, 3, 1)
        self.singleProcessDate.setText(self.date)
        self.singleProcessDate.setEnabled(False)

        # Label for manual single cluster
        self.singleProcessClusterLabel = QtGui.QLabel(self)
        self.singleProcessClusterLabel.resize(
            self.singleProcessClusterLabel.sizeHint())
        self.singleProcessClusterLabel.setText("Cluster")
        self.mainGrid.addWidget(self.singleProcessClusterLabel, 4, 0)

        # Input for manual single cluster
        self.singleProcessCluster = QtGui.QLineEdit(self)
        self.singleProcessCluster.setPlaceholderText("Name")
        self.singleProcessCluster.textChanged.connect(
            self.SingleProcessClusterChange)
        self.mainGrid.addWidget(self.singleProcessCluster, 4, 1)
        self.singleProcessCluster.setText(self.cluster)

        # Label for coordinate tolerance
        self.cooTolInputLabel = QtGui.QLabel(self)
        self.cooTolInputLabel.resize(self.cooTolInputLabel.sizeHint())
        self.cooTolInputLabel.setText("Coo. Tol.")
        self.mainGrid.addWidget(self.cooTolInputLabel, 5, 0)

        # Input for coordinate tolerance
        self.cooTolInput = QtGui.QLineEdit(self)
        self.cooTolInput.setPlaceholderText("Tolerance")
        self.cooTolInput.setText(str(self.cooTol))
        self.cooTolInput.resize(self.cooTolInput.sizeHint())
        self.cooTolInput.textChanged.connect(self.CooTolChange)
        self.mainGrid.addWidget(self.cooTolInput, 5, 1)

        # Label for magnitude tolerance
        self.magTolInputLabel = QtGui.QLabel(self)
        self.magTolInputLabel.resize(self.magTolInputLabel.sizeHint())
        self.magTolInputLabel.setText("Mag. Tol.")
        self.mainGrid.addWidget(self.magTolInputLabel, 6, 0)

        # Input for magnitude tolerance
        self.magTolInput = QtGui.QLineEdit(self)
        self.magTolInput.setPlaceholderText("Tolerance")
        self.magTolInput.setText(str(self.magTol))
        self.magTolInput.resize(self.magTolInput.sizeHint())
        self.magTolInput.textChanged.connect(self.MagTolChange)
        self.mainGrid.addWidget(self.magTolInput, 6, 1)

        # Label for manual constant threshold
        self.thresholdInputLabel = QtGui.QLabel(self)
        self.thresholdInputLabel.resize(self.thresholdInputLabel.sizeHint())
        self.thresholdInputLabel.setText("R-H Threshold")
        self.mainGrid.addWidget(self.thresholdInputLabel, 7, 0)
        self.thresholdInputLabel.setEnabled(False)

        # Input for manual constant threshold
        self.thresholdInput = QtGui.QLineEdit(self)
        self.thresholdInput.setPlaceholderText("Threshold")
        self.thresholdInput.setText(str(self.threshold))
        self.thresholdInput.resize(self.thresholdInput.sizeHint())
        self.thresholdInput.textChanged.connect(self.ThresholdInputChange)
        self.mainGrid.addWidget(self.thresholdInput, 7, 1)
        self.thresholdInput.setEnabled(False)

        # Label for B-V maximum
        self.B_VMaxInputLabel = QtGui.QLabel(self)
        self.B_VMaxInputLabel.resize(self.B_VMaxInputLabel.sizeHint())
        self.B_VMaxInputLabel.setText("B-V Max")
        self.mainGrid.addWidget(self.B_VMaxInputLabel, 8, 0)

        # Input for B-V maximum
        self.B_VMaxInput = QtGui.QLineEdit(self)
        self.B_VMaxInput.setPlaceholderText("Max")
        self.B_VMaxInput.setText('%.3f' % self.B_VMax)
        self.B_VMaxInput.resize(self.B_VMaxInput.sizeHint())
        self.B_VMaxInput.textChanged.connect(self.B_VMaxInputChange)
        self.mainGrid.addWidget(self.B_VMaxInput, 8, 1)

        # Label for B-V minimum
        self.B_VMinInputLabel = QtGui.QLabel(self)
        self.B_VMinInputLabel.resize(self.B_VMinInputLabel.sizeHint())
        self.B_VMinInputLabel.setText("B-V Min")
        self.mainGrid.addWidget(self.B_VMinInputLabel, 9, 0)

        # Input for B-V minimum
        self.B_VMinInput = QtGui.QLineEdit(self)
        self.B_VMinInput.setPlaceholderText("Min")
        self.B_VMinInput.setText('%.3f' % self.B_VMin)
        self.B_VMinInput.resize(self.B_VMinInput.sizeHint())
        self.B_VMinInput.textChanged.connect(self.B_VMinInputChange)
        self.mainGrid.addWidget(self.B_VMinInput, 9, 1)

        # Process button
        self.processButton = QtGui.QPushButton("Process", self)
        self.processButton.resize(self.processButton.sizeHint())
        self.processButton.clicked.connect(self.Process)
        self.mainGrid.addWidget(self.processButton, 9, 2, 1, 3)

        # Console output
        self.logOutput = QtGui.QTextEdit(self)
        self.logOutput.setReadOnly(True)
        self.mainGrid.addWidget(self.logOutput, 2, 2, 7, 3)

        self.show()

    def ProcessTypeChange(self, text):
        """Controls GUI and class elements when processType changes."""
        self.process_type = text

        if text == "Single Date":
            self.singleProcessDate.setEnabled(True)
            self.singleProcessDateLabel.setEnabled(True)
            self.singleProcessCluster.setEnabled(True)
            self.singleProcessClusterLabel.setEnabled(True)
            self.singleProcessDate.setText(self.date)
            self.singleProcessCluster.setText(self.cluster)
        elif text == "Single Cluster":
            self.singleProcessDate.setEnabled(False)
            self.singleProcessDateLabel.setEnabled(False)
            self.singleProcessCluster.setEnabled(True)
            self.singleProcessClusterLabel.setEnabled(True)
            self.singleProcessCluster.setText(self.cluster)
        elif text == "Full":
            self.singleProcessDate.setEnabled(False)
            self.singleProcessDateLabel.setEnabled(False)
            self.singleProcessCluster.setEnabled(False)
            self.singleProcessClusterLabel.setEnabled(False)

    def ThresholdTypeChange(self, text):
        """Controls GUI and class elements when threshold type changes."""
        self.threshold_type = text

    def SingleProcessDateChange(self, text):
        """Controls GUI and class elements when the manual process date
           changes."""
        self.date = text

    def SingleProcessClusterChange(self, text):
        """Controls GUI and class elements when the manual process cluster
           changes."""
        self.cluster = text

    def AutoThresholdCheckChange(self, state):
        """Controls GUI and class elements when the auto-threshold checkbox
           changes."""
        if state == QtCore.Qt.Checked:
            self.thresholdInputLabel.setEnabled(False)
            self.thresholdInput.setEnabled(False)
        else:
            self.thresholdInputLabel.setEnabled(True)
            self.thresholdInput.setEnabled(True)

    def CooTolChange(self, text):
        """Controls GUI and class elements when the coordinate tolerance input
           changes."""
        try:
            self.cooTol = float(text)
        except Exception:
            pass

    def MagTolChange(self, text):
        """Controls GUI and class elements when the magnitude tolerance input
           changes."""
        try:
            self.magTol = float(text)
        except Exception:
            pass

    def ThresholdInputChange(self, text):
        """Controls GUI and class elements when the manual threshold input
           changes."""
        try:
            self.threshold = float(text)
        except Exception:
            pass

    def B_VMaxInputChange(self, text):
        """Controls GUI and class elements when the B-V maximum value input
           changes."""
        try:
            self.B_VMax = float(text)
        except Exception:
            pass

    def B_VMinInputChange(self, text):
        """Controls GUI and class elements when the B-V minimum value input
           changes."""
        try:
            self.B_VMin = float(text)
        except Exception:
            pass

    def Process(self):
        """Calls the process controller when the 'Process' button is called."""
        Process(self)
