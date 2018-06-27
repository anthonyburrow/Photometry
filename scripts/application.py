from PyQt4 import QtGui, QtCore
from process_control import ProcessControl


class Application(QtGui.QMainWindow):

    def __init__(self):
        super(Application, self).__init__()

        # Set defaults
        self.process_type = "Single"
        self.phot_type = "psf"
        self.threshold_type = "Linear"
        self.date = None
        self.cluster = None
        self.cooTol = 5
        self.magTol = 0.5
        self.threshold = -3.75
        self.B_VMax = 1.00
        self.B_VMin = -0.1

        # Configure window
        # self.setGeometry(50, 50, 600, 400)
        self.setWindowTitle("Photometry")
        # self.setWindowIcon(QtGui.QIcon('logo'))

        # Configure grid layouts
        centralWidget = QtGui.QWidget()
        self.mainGrid = QtGui.QGridLayout()
        centralWidget.setLayout(self.mainGrid)
        self.setCentralWidget(centralWidget)

        # Photometry type
        self.photType = QtGui.QComboBox(self)
        self.photType.addItem("PSF + Aperture")
        self.photType.addItem("Aperture")
        self.photType.resize(self.photType.sizeHint())
        self.photType.activated[str].connect(self.PhotTypeChange)
        self.mainGrid.addWidget(self.photType, 0, 0)

        # Process type
        self.processType = QtGui.QComboBox(self)
        self.processType.addItem("Single")
        self.processType.addItem("Full")
        self.processType.resize(self.processType.sizeHint())
        self.processType.activated[str].connect(self.ProcessTypeChange)
        self.mainGrid.addWidget(self.processType, 0, 1)

        # Threshold type
        self.thresholdType = QtGui.QComboBox(self)
        self.thresholdType.addItem("Linear")
        self.thresholdType.addItem("Constant")
        self.thresholdType.resize(self.thresholdType.sizeHint())
        self.thresholdType.activated[str].connect(self.ThresholdTypeChange)
        self.mainGrid.addWidget(self.thresholdType, 0, 2)

        # Auto threshold checkbox
        self.autoThresholdCheck = QtGui.QCheckBox("Auto Threshold", self)
        self.autoThresholdCheck.resize(self.autoThresholdCheck.sizeHint())
        self.autoThresholdCheck.toggle()
        self.autoThresholdCheck.stateChanged.connect(self.AutoThresholdCheckChange)
        self.mainGrid.addWidget(self.autoThresholdCheck, 0, 3)

        # Run match process
        self.matchCheck = QtGui.QCheckBox("Match", self)
        self.matchCheck.resize(self.matchCheck.sizeHint())
        self.mainGrid.addWidget(self.matchCheck, 1, 0)

        # Run Be filter process
        self.befilterCheck = QtGui.QCheckBox("Filter Be Candidates", self)
        self.befilterCheck.resize(self.befilterCheck.sizeHint())
        self.mainGrid.addWidget(self.befilterCheck, 1, 1)

        # Run plot process
        self.plotCheck = QtGui.QCheckBox("Plot", self)
        self.plotCheck.resize(self.plotCheck.sizeHint())
        self.plotCheck.stateChanged.connect(self.PlotCheckChange)
        self.mainGrid.addWidget(self.plotCheck, 1, 2)

        # Plot CMDs
        self.plotCMDCheck = QtGui.QCheckBox("CMD", self)
        self.plotCMDCheck.resize(self.plotCMDCheck.sizeHint())
        self.plotCMDCheck.setEnabled(False)
        self.mainGrid.addWidget(self.plotCMDCheck, 1, 3)

        # Plot 2CDs
        self.plot2CDCheck = QtGui.QCheckBox("2CD", self)
        self.plot2CDCheck.resize(self.plot2CDCheck.sizeHint())
        self.plot2CDCheck.setEnabled(False)
        self.mainGrid.addWidget(self.plot2CDCheck, 1, 4)

        # Run scale process
        self.scaleCheck = QtGui.QCheckBox("Scale", self)
        self.scaleCheck.resize(self.scaleCheck.sizeHint())
        self.mainGrid.addWidget(self.scaleCheck, 2, 0)

        # Label for manual single date
        self.singleProcessDateLabel = QtGui.QLabel(self)
        self.singleProcessDateLabel.resize(self.singleProcessDateLabel.sizeHint())
        self.singleProcessDateLabel.setText("Date")
        self.mainGrid.addWidget(self.singleProcessDateLabel, 3, 0)

        # Input for manual single date
        self.singleProcessDate = QtGui.QLineEdit(self)
        self.singleProcessDate.setMaxLength(8)
        self.singleProcessDate.setPlaceholderText("MMDDYYYY")
        self.singleProcessDate.textChanged.connect(self.SingleProcessDateChange)
        self.mainGrid.addWidget(self.singleProcessDate, 3, 1)

        # Label for manual single cluster
        self.singleProcessClusterLabel = QtGui.QLabel(self)
        self.singleProcessClusterLabel.resize(self.singleProcessClusterLabel.sizeHint())
        self.singleProcessClusterLabel.setText("Cluster")
        self.mainGrid.addWidget(self.singleProcessClusterLabel, 4, 0)

        # Input for manual single cluster
        self.singleProcessCluster = QtGui.QLineEdit(self)
        self.singleProcessCluster.setPlaceholderText("Name")
        self.singleProcessCluster.textChanged.connect(self.SingleProcessClusterChange)
        self.mainGrid.addWidget(self.singleProcessCluster, 4, 1)

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
        self.B_VMaxInputLabel.setText("B-V Min")
        self.mainGrid.addWidget(self.B_VMaxInputLabel, 8, 0)

        # Input for B-V maximum
        self.B_VMaxInput = QtGui.QLineEdit(self)
        self.B_VMaxInput.setPlaceholderText("Max")
        self.B_VMaxInput.setText(str(self.B_VMax))
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
        self.B_VMinInput.setText(str(self.B_VMin))
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
        """Controls GUI and class elements when phot type changes."""
        if text == "PSF + Aperture":
            self.phot_type = "psf"
        elif text == "Aperture":
            self.phot_type = "aperture"

    def ThresholdTypeChange(self, text):
        """Controls GUI and class elements when threshold type changes."""
        self.threshold_type = text

    def SingleProcessDateChange(self, text):
        """Controls GUI and class elements when the manual process date changes."""
        self.date = text

    def SingleProcessClusterChange(self, text):
        """Controls GUI and class elements when the manual process cluster changes."""
        self.cluster = text

    def AutoThresholdCheckChange(self, state):
        """Controls GUI and class elements when the auto-threshold checkbox changes."""
        if state == QtCore.Qt.Checked:
            self.thresholdInputLabel.setEnabled(False)
            self.thresholdInput.setEnabled(False)
        else:
            self.thresholdInputLabel.setEnabled(True)
            self.thresholdInput.setEnabled(True)

    def CooTolChange(self, text):
        """Controls GUI and class elements when the coordinate tolerance input changes."""
        try:
            self.cooTol = float(text)
        except Exception:
            pass

    def MagTolChange(self, text):
        """Controls GUI and class elements when the magnitude tolerance input changes."""
        try:
            self.magTol = float(text)
        except Exception:
            pass

    def ThresholdInputChange(self, text):
        """Controls GUI and class elements when the manual threshold input changes."""
        try:
            self.threshold = float(text)
        except Exception:
            pass

    def B_VMaxInputChange(self, text):
        """Controls GUI and class elements when the B-V maximum value input changes."""
        try:
            self.B_VMax = float(text)
        except Exception:
            pass

    def B_VMinInputChange(self, text):
        """Controls GUI and class elements when the B-V minimum value input changes."""
        try:
            self.B_VMin = float(text)
        except Exception:
            pass

    def PlotCheckChange(self, state):
        """Controls GUI and class elements when the plot checkbox changes."""
        if state == QtCore.Qt.Checked:
            self.plotCMDCheck.setEnabled(True)
            self.plot2CDCheck.setEnabled(True)
        else:
            self.plotCMDCheck.setEnabled(False)
            self.plot2CDCheck.setEnabled(False)

    def Process(self):
        """Calls the process controller when the 'Process' button is called."""
        processControl = ProcessControl(self)
        processControl.Process()
