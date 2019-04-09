#!/usr/bin/python
# -*- coding: utf-8 -*-

"""
ZetCode PyQt4 tutorial 

This program creates a menubar. The
menubar has one menu with an exit action.

author: Jan Bodnar
website: zetcode.com 
last edited: August 2011
"""

import sys
from PyQt4 import QtGui, QtCore
from util import Time

class Example(QtGui.QMainWindow):
    
    def __init__(self, parent=None):
        super(Example, self).__init__(parent)
        
        newAction = QtGui.QAction('&New', self)        
        newAction.setShortcut('Ctrl+N')
        newAction.setStatusTip('Create new file')
        newAction.triggered.connect(self.new)
        saveAction = QtGui.QAction('&Save', self)        
        saveAction.setShortcut('Ctrl+S')
        saveAction.setStatusTip('Save current file')
        saveAction.triggered.connect(self.save)
        openAction = QtGui.QAction('&Open', self)        
        openAction.setShortcut('Ctrl+O')
        openAction.setStatusTip('Open an existing file')
        openAction.triggered.connect(self.open)
        exitAction = QtGui.QAction('&Exit', self)        
        exitAction.setShortcut('Ctrl+E')
        exitAction.setStatusTip('Exit the schedule')
        exitAction.triggered.connect(QtGui.qApp.quit)

        self.statusBar()

        menubar = self.menuBar()
        fileMenu = menubar.addMenu('&File')
        fileMenu.addAction(newAction)
        fileMenu.addAction(openAction)
        fileMenu.addAction(saveAction)
        fileMenu.addAction(exitAction)
        self.form_widget = FormWidget(self)
        self.setCentralWidget(self.form_widget)
        self.show()
        self.timer = QtCore.QTimer(self)
        
        self.timer.singleShot(1,self.inc_time)

    def new(self):
        print 'New File was selected'
        
    def open(self):
        print 'Open File was selected'
        
    def save(self):
        print 'Save File was selected'
        
    def inc_time(self):
        self.timer.singleShot(1000,self.inc_time)
        t = Time.now()
        self.form_widget.label.setText(t.iso[:19])
        
class FormWidget(QtGui.QWidget):
    def __init__(self, parent):
        super(FormWidget,self).__init__(parent)
        # H-box for Time widget
        self.label = QtGui.QLabel("") #bg="yellow",font="Helvetica 16 bold"
        self.label.setStyleSheet('font-size: 16pt; font-weight: bold; background-color: yellow')
        timeframe = QtGui.QHBoxLayout()
        timeframe.addStretch(1)
        timeframe.addWidget(self.label)
        timeframe.addStretch(1)
        
        self.error = ''   # Error string to be added to source on front panel

        layout = QtGui.QHBoxLayout()
        self.b1 = QtGui.QRadioButton("D")
        self.b1.setChecked(True)
        self.b_text = "D"
        self.b1.toggled.connect(lambda:self.btnstate(self.b1))
        layout.addWidget(self.b1)
        self.b2 = QtGui.QRadioButton("H")
        self.b2.toggled.connect(lambda:self.btnstate(self.b2))
        layout.addWidget(self.b2)
        self.b3 = QtGui.QRadioButton("M")
        self.b3.toggled.connect(lambda:self.btnstate(self.b3))
        layout.addWidget(self.b3)
        self.b4 = QtGui.QRadioButton("S")
        self.b4.toggled.connect(lambda:self.btnstate(self.b4))
        layout.addWidget(self.b4)
        self.source = QtGui.QLabel("            Source: None        "+"Phase Tracking: False")
        layout.addWidget(self.source)
        layout.addStretch(1)
        vbox = QtGui.QVBoxLayout()
        vbox.addLayout(timeframe)
        vbox.addLayout(layout)
        vbox.addStretch(1)
        
        self.setLayout(vbox)
        self.setGeometry(300, 300, 300, 200)
        self.setWindowTitle("RadioButton demo")
        
    def btnstate(self,b):
        if b.isChecked():
            self.source.setText("            Source: "+b.text()+"      "+"Phase Tracking: False")
            #print b.text()+' is selected'
            self.b_text = b.text()

def main():
    
    app = QtGui.QApplication(sys.argv)
    ex = Example()
    sys.exit(app.exec_())


if __name__ == '__main__':
    main()    