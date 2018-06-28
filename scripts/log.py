import sys
from PyQt4 import QtGui


class Log(object):
    def __init__(self, edit, logfile):
        self.out = sys.stdout
        self.textEdit = edit
        self.logfile = logfile

    def write(self, message):
        self.out.write(message)
        self.logfile.write(message)
        self.textEdit.insertPlainText(message)
        self.textEdit.moveCursor(QtGui.QTextCursor.End)

    def flush(self):
        self.out.flush()
