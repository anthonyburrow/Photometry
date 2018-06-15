from application import Application
import sys
from PyQt4 import QtGui
from make_files import MakeFiles

root = "../photometry/"

makeFiles = MakeFiles(root)
makeFiles.ObsList()

print ("Opening application...")

app = QtGui.QApplication(sys.argv)
gui = Application()
sys.exit(app.exec_())