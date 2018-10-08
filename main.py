import sys
import os.path

from PyQt4 import QtGui

from Analysis.application import Application
from Analysis.log import Log

print("Opening application...")

app = QtGui.QApplication(sys.argv)
gui = Application()

# Set up directories for processing
setup_directories = [
    'output/',
    'photometry/'
]
for directory in setup_directories:
    if not os.path.exists(directory):
        os.makedirs(directory)

# Log and output
F = open('output/log_last_run.txt', 'w')
sys.stdout = Log(gui.logOutput, F)

sys.exit(app.exec_())
F.close()
