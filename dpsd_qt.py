#!/usr/bin/env python

__author__  = 'Giovanni Tardini (Tel. 1898)'
__version__ = '0.0.1'
__date__    = '29.03.2022'

import os, sys, logging

from PyQt5 import uic
from PyQt5.QtWidgets import QMainWindow, QWidget, QApplication, QVBoxLayout, QHBoxLayout, QGridLayout, QMenu, QAction, QLabel, QPushButton, QLineEdit, QRadioButton, QCheckBox, QButtonGroup
from PyQt5.QtGui import QPixmap, QIcon
from PyQt5.QtCore import Qt, QRect, QSize

import numpy as np
import matplotlib.pylab as plt
import dixm, dpsd, plot_dpsd
try:
    import aug_sfutils as sf
except:
    pass


fmt = logging.Formatter('%(asctime)s | %(name)s | %(levelname)s: %(message)s', '%H:%M:%S')
hnd = logging.StreamHandler()
hnd.setFormatter(fmt)
logger = logging.getLogger('DPSD_GUI')
logger.addHandler(hnd)
logger.setLevel(logging.DEBUG)

frc = '#b0d0b0'  # Frame, Notebook, Checkbutton
tbc = '#eaeaea'  # Toolbar, Button
figc = '#70a0c0' # Figure
colors  = 2*['#c00000', '#00c000', '#0000c0', '#b0b000', '#b000b0', '#00b0b0']

dpsd_dir = os.path.dirname(os.path.realpath(__file__))


class DPSD(QMainWindow):


    def __init__(self):

        super().__init__()
        self.setWindowTitle('DPSD')

        xwin  = 640
        ywin  = 630
        yhead = 44
        ybar  = 50

        head = QWidget(self)
        tbar = QWidget(self)
        body = QWidget(self)
        head.setGeometry(QRect(0,     0, xwin, yhead))
        tbar.setGeometry(QRect(0, yhead, xwin, ybar))
        body.setGeometry(QRect(0, yhead+ybar, xwin, ywin-yhead-ybar))
        header_grid = QGridLayout(head) 
        tbar_grid   = QGridLayout(tbar) 
        entry_grid  = QGridLayout(body)

# Menubar

        menubar = self.menuBar()
        fileMenu = QMenu('&File', self)
        helpMenu = QMenu('&Help', self)
        menubar.addMenu(fileMenu)
        menubar.addMenu(helpMenu)

        runAction  = QAction('&Run'     , fileMenu)
        plotAction = QAction('&Plot'    , fileMenu)
        wsfAction  = QAction('&Write SF', fileMenu)
        exitAction = QAction('&Exit'    , fileMenu)
        runAction.triggered.connect(self.run)
        plotAction.triggered.connect(self.plot)
        wsfAction.triggered.connect(self.write_sf)
        exitAction.triggered.connect(sys.exit)
        fileMenu.addAction(runAction)
        fileMenu.addAction(plotAction)
        fileMenu.addAction(wsfAction)
        fileMenu.addSeparator()
        fileMenu.addAction(exitAction)

        aboutAction = QAction('&About', helpMenu)
        aboutAction.triggered.connect(self.about)
        helpMenu.addAction(aboutAction)

        header_grid.addWidget(menubar, 0, 0, 1, 10)

# Icons
        dum_lbl  = QLabel(600*' ')
        fmap = {'play': self.run, 'plot': self.plot, \
            'save': self.write_sf, 'exit': sys.exit}
        for jpos, lbl in enumerate(fmap.keys()):
            but = QPushButton()
            but.setIcon(QIcon('%s/%s.gif' %(dpsd_dir, lbl)))
            but.setIconSize(QSize(ybar, ybar))
            but.clicked.connect(fmap[lbl])
            tbar_grid.addWidget(but, 0, jpos)
        tbar_grid.addWidget(dum_lbl,  0, 4, 1, 10)

# User options

        setup_en_d = dixm.DIX().xml2dict('/afs/ipp/home/g/git/DPSD/xml/default.xml')
        self.gui = {}
        user = os.getenv('USER')

# Entry widgets
        jcol = 0
        jrow = 0
        key = 'Path'
        self.gui[key] = QLineEdit(setup_en_d[key])
        entry_grid.addWidget(QLabel(key), jrow, 0)
        entry_grid.addWidget(self.gui[key], jrow, 1, 1, 3)
        jrow += 1
        for key, val in setup_en_d.items():
            if key == 'Path':
                continue
            label = QLabel(key)
            self.gui[key] = QLineEdit(val)
            self.gui[key].setFixedWidth(90)
            if jrow == 18:
                jrow = 1
                jcol += 2
            entry_grid.addWidget(label        , jrow, jcol)
            entry_grid.addWidget(self.gui[key], jrow, jcol+1)
            jrow += 1

# Checkbutton

        keys = ['SubtBaseline', 'LEDcorrection']
        cb_d = {'SubtBaseline': 'Subtract baseline', 'LEDcorrection': 'LED correction'}

        jrow = 18
        for key, lbl in cb_d.items():
            jrow += 1
            self.gui[key] = QCheckBox(lbl)
            entry_grid.addWidget(self.gui[key], jrow, 0, 1, 2)
            self.gui[key].setChecked(True)

        self.setStyleSheet("QLabel { width: 4 }")
        self.setStyleSheet("QLineEdit { width: 4 }")
        self.setGeometry(10, 10, xwin, ywin)
        self.setWindowTitle('Confinement')
        self.show()


    def about(self):

        mytext = 'Documentation at <a href="http://www.aug.ipp.mpg.de/~git/tot/index.html">TOT/TTH diagnostic homepage</a>'
        h = tkhyper.HyperlinkMessageBox("Help", mytext, "500x60")


    def run(self):

        dpsd_dic = {}
        for key, val in self.gui.items():
            if isinstance(val, QLineEdit):
                dpsd_dic[key] = val.text()
            elif isinstance(val, QCheckBox):
                dpsd_dic[key] = val.isChecked()
        self.dp = dpsd.DPSD(dpsd_dic)


    def plot(self):

        if not hasattr(self, 'dp'):
            logger.error('Run code before plotting')
            return

        self.wid = plot_dpsd.plotWindow()
        fig1 = plot_dpsd.fig_pha(self.dp)
        fig2 = plot_dpsd.fig_phs(self.dp)
        fig3 = plot_dpsd.fig_cnt(self.dp)
        fig4 = plot_dpsd.fig_pmg(self.dp)

        self.wid.addPlot('PH-PS separation', fig1)
        self.wid.addPlot('PHA spectrum'    , fig2)
        self.wid.addPlot('Count rates'     , fig3)
        self.wid.addPlot('PM gain'         , fig4)

        self.wid.show()
 

    def write_sf(self):

        if hasattr(self, 'dp'):
            self.dp.sfwrite(exp='git')
        else:
            logger.error('Run DPSD first, then write Shotfile')


if __name__ == '__main__':


    app = QApplication(sys.argv)
    main = DPSD()
    app.exec()
