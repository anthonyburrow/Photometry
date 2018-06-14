import sys
from PyQt4 import QtGui, QtCore

class Application(QtGui.QMainWindow):

	def __init__(self):
		super(Application, self).__init__()

		self.setGeometry(50, 50, 500, 300)
		self.setWindowTitle("Photometry")
		#self.setWindowIcon(QtGui.QIcon('logo'))

		self.home()

	def home(self):
		btn = QtGui.QPushButton("Quit", self)
		btn.clicked.connect(QtCore.QCoreApplication.instance().quit)

		btn.resize(100, 100)
		btn.move(100, 100)

		self.show()