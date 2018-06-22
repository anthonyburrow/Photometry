from application import Application
import sys
from PyQt4 import QtGui
from log import Log

print("Opening application...")

app = QtGui.QApplication(sys.argv)
gui = Application()

sys.stdout = Log(gui.logOutput)

sys.exit(app.exec_())
