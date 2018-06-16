import sys
from PyQt4 import QtGui, QtCore


class Log(object):
    def __init__(self, edit):
        self.out = sys.stdout
        self.textEdit = edit

    def write(self, message):
        self.out.write(message)
        self.textEdit.insertPlainText(message)
        self.textEdit.moveCursor(QtGui.QTextCursor.End)

    def flush(self):
        self.out.flush()
