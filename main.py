import sys

from PyQt4 import QtGui

from Analysis.application import Application
from Analysis.log import Log
from Analysis.initialize import Setup


def main():
    print("Opening application...")

    app = QtGui.QApplication(sys.argv)
    gui = Application()

    Setup()

    # Log and output
    F = open('output/log_last_run.txt', 'w')
    sys.stdout = Log(gui.logOutput, F)

    sys.exit(app.exec_())
    F.close()


if __name__ == '__main__':
    main()
