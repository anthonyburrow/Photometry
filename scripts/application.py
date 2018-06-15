from PyQt4 import QtGui, QtCore
from match import Match
from be_filter import BeFilter
from plot import Plot
from scale import Scale

class Application(QtGui.QMainWindow):

	def __init__(self):
		super(Application, self).__init__()

		self.setGeometry(50, 50, 500, 300)
		self.setWindowTitle("Photometry")
		#self.setWindowIcon(QtGui.QIcon('logo'))

		self.process_type = "Single"
		self.phot_type = "psf"
		self.date = None
		self.cluster = None

		self.processType = QtGui.QComboBox(self)
		self.photType = QtGui.QComboBox(self)
		self.singleProcessDate = QtGui.QLineEdit(self)
		self.matchCheck = QtGui.QCheckBox("Match", self)
		self.befilterCheck = QtGui.QCheckBox("Filter Be Candidates", self)
		self.plotCheck = QtGui.QCheckBox("Plot", self)
		self.processButton = QtGui.QPushButton("PROCESS", self)

		self.home()

	def home(self):
		self.processType.move(0, 0)
		self.processType.addItem("Single")
		self.processType.addItem("Full")
		self.processType.activated[str].connect(self.ProcessTypeChange)

		self.photType.move(0, 100)
		self.photType.addItem("PSF + Aperture")
		self.photType.addItem("Aperture")
		self.photType.activated[str].connect(self.PhotTypeChange)

		self.singleProcessDate.move(150, 0)
		self.singleProcessDate.setMaxLength(8)
		self.singleProcessDate.setPlaceholderText("MMDDYYYY")
		self.singleProcessDate.textChanged.connect(self.SingleProcessDateChange)

		self.matchCheck.move(200, 200)
		self.matchCheck.resize(self.matchCheck.sizeHint())
		self.befilterCheck.move(250, 200)
		self.befilterCheck.resize(self.befilterCheck.sizeHint())
		self.plotCheck.move(367, 200)
		self.plotCheck.resize(self.plotCheck.sizeHint())

		self.processButton.move(200, 300)
		self.processButton.resize(self.processButton.sizeHint())
		self.processButton.clicked.connect(self.ProcessControl)

		self.show()

	def ProcessTypeChange(self, text):
		self.process_type = text

		if text == "Single":
			self.singleProcessDate.setEnabled(True)
			self.singleProcessDate.setText(self.date)
		elif text == "Full":
			self.singleProcessDate.setEnabled(False)

	def PhotTypeChange(self, text):
		if text == "PSF + Aperture":
			self.phot_type = "psf"
		elif text == "Aperture":
			self.phot_type = "aperture"

	def SingleProcessDateChange(self, text):
		self.date = text

	def ProcessControl(self):
		option = str(processType.currentText())

		if option == "Single":
			ProcessDate(self.cluster, self.date)

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
			print ("Compiling all data for " + self.cluster + " on " + self.date + "...")
			match = Match(output_directory, input_directory, self.phot_type, 5.0, 0.5)
			match.LowError()

		if self.matchCheck.isChecked():
			print ("Extracting Be candidate data...")
			beFilter = BeFilter(output_directory, self.phot_type, true)
			beFilter.LowError()

		if self.matchCheck.isChecked():
			print ("Generating plots...")
			plot = Plot(output_directory, true, true)
			plot.ColorMagnitudeDiagram()
			plot.TwoColorDiagram()

	def ProcessCluster(self):

		print("Scaling observation nights for " + self.cluster)
		scale = Scale(self.cluster, self.phot_type, 10)