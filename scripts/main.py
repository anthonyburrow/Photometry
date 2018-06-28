from application import Application
import sys
from PyQt4 import QtGui
from log import Log
import os.path

print("Opening application...")

app = QtGui.QApplication(sys.argv)
gui = Application()

# Log and output
if not os.path.exists("../output/"):
    os.makedirs("../output/")

F = open('../output/log_last_run.txt', 'w')
sys.stdout = Log(gui.logOutput, F)

sys.exit(app.exec_())
F.close()
