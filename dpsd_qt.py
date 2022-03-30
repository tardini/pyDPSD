#!/usr/bin/env python

__author__  = 'Giovanni Tardini (Tel. 1898)'
__version__ = '0.0.1'
__date__    = '29.03.2022'

import os, sys, logging

try:
    from PyQt5.QtWidgets import QMainWindow, QWidget, QApplication, QGridLayout, QMenu, QAction, QLabel, QPushButton, QLineEdit, QCheckBox, QFileDialog
    from PyQt5.QtGui import QPixmap, QIcon
    from PyQt5.QtCore import Qt, QRect, QSize
    qt5 = True
except:
    from PyQt4.QtCore import Qt, QRect, QSize
    from PyQt4.QtGui import QPixmap, QIcon, QMainWindow, QWidget, QApplication, QGridLayout, QMenu, QAction, QLabel, QPushButton, QLineEdit, QCheckBox, QFileDialog
    qt5 = False

import numpy as np
import dixm, dpsd, plot_dpsd
try:
    import aug_sfutils as sf
except:
    pass

xml = dixm.DIX()

fmt = logging.Formatter('%(asctime)s | %(name)s | %(levelname)s: %(message)s', '%H:%M:%S')
hnd = logging.StreamHandler()
hnd.setFormatter(fmt)
logger = logging.getLogger('DPSD_GUI')
logger.addHandler(hnd)
logger.setLevel(logging.INFO)

frc = '#b0d0b0'  # Frame, Notebook, Checkbutton
tbc = '#eaeaea'  # Toolbar, Button
figc = '#70a0c0' # Figure
colors  = 2*['#c00000', '#00c000', '#0000c0', '#b0b000', '#b000b0', '#00b0b0']

dpsd_dir = os.path.dirname(os.path.realpath(__file__))


class DPSD(QMainWindow):


    def __init__(self):

        if qt5:
            super().__init__()
        else:
            super(QMainWindow, self).__init__()

        self.setWindowTitle('DPSD')

        xwin  = 450
        ywin  = 660
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
        xmlMenu  = QMenu('&Setup', self)
        helpMenu = QMenu('&Help', self)
        menubar.addMenu(fileMenu)
        menubar.addMenu(xmlMenu)
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

        loadAction  = QAction('&Load Setup', xmlMenu)
        saveAction  = QAction('&Save Setup', xmlMenu)
        xmlMenu.addAction(loadAction)
        xmlMenu.addAction(saveAction)
        loadAction.triggered.connect(self.load_xml)
        saveAction.triggered.connect(self.save_xml)

        aboutAction = QAction('&About', helpMenu)
        aboutAction.triggered.connect(self.about)
        helpMenu.addAction(aboutAction)

        header_grid.addWidget(menubar, 0, 0, 1, 10)

# Icons
        dum_lbl  = QLabel(600*' ')
        fmap = {'play': self.run, 'plot': self.plot, \
            'save': self.write_sf, 'exit': sys.exit}
        for jpos, lbl in enumerate(['play', 'exit', 'plot', 'save']):
            but = QPushButton()
            but.setIcon(QIcon('%s/%s.gif' %(dpsd_dir, lbl)))
            but.setIconSize(QSize(ybar, ybar))
            but.clicked.connect(fmap[lbl])
            tbar_grid.addWidget(but, 0, jpos)
        tbar_grid.addWidget(dum_lbl,  0, 4, 1, 10)

# User options

        setup_en_d = xml.xml2dict('%s/xml/default.xml' %dpsd_dir)
        self.gui = {}
        user = os.getenv('USER')

# Entry widgets
        cb_d = {'SubtBaseline': 'Subtract baseline', 'LEDcorrection': 'LED correction', 'SFwrite': 'Write shotfiles'}
        jcol = 0
        jrow = 0
        key = 'HAfile'
        self.gui[key] = QLineEdit(setup_en_d[key])
        entry_grid.addWidget(QLabel(key), jrow, 0)
        entry_grid.addWidget(self.gui[key], jrow, 1, 1, 3)
        key = 'Shots'
        self.gui[key] = QLineEdit(setup_en_d[key])
        entry_grid.addWidget(QLabel(key), jrow+1, 0)
        entry_grid.addWidget(self.gui[key], jrow+1, 1, 1, 3)
        row_init = 2
        row_end = 20
        jrow += row_init
        for key, val in setup_en_d.items():
            if key in ('HAfile', 'Shots'):
                continue
            if key in cb_d.keys():
                continue
            label = QLabel(key)
            self.gui[key] = QLineEdit(val)
            self.gui[key].setFixedWidth(90)
            if jrow == row_end:
                jrow = row_init
                jcol += 2
            entry_grid.addWidget(label        , jrow, jcol)
            entry_grid.addWidget(self.gui[key], jrow, jcol+1)
            jrow += 1

# Checkbutton

        jrow = row_end
        for key, lbl in cb_d.items():
            jrow += 1
            self.gui[key] = QCheckBox(lbl)
            entry_grid.addWidget(self.gui[key], jrow, 0, 1, 2)
            if setup_en_d[key].lower().strip() == 'true':
                self.gui[key].setChecked(True)

        self.setStyleSheet("QLabel { width: 4 }")
        self.setStyleSheet("QLineEdit { width: 4 }")
        self.setGeometry(10, 10, xwin, ywin)
        self.setWindowTitle('DPSD')
        self.show()


    def about(self):

        mytext = 'Documentation at <a href="http://www.aug.ipp.mpg.de/~git/tot/index.html">TOT/TTH diagnostic homepage</a>'
        h = tkhyper.HyperlinkMessageBox("Help", mytext, "500x60")


    def get_gui(self):

        dpsd_dic = {}
        for key, val in self.gui.items():
            if isinstance(val, QLineEdit):
                dpsd_dic[key] = val.text()
            elif isinstance(val, QCheckBox):
                dpsd_dic[key] = val.isChecked()

        return dpsd_dic


    def set_gui(self, xml_d):

        for key, val in xml_d.items():
            val = val.strip()
            val_low = val.lower()
            if val_low == 'false':
                self.gui[key].setChecked(False)
            elif val_low == 'true':
                self.gui[key].setChecked(True)
            elif val_low == '':
                self.gui[key].setText(' ')
            else:
                self.gui[key].setText(val)


    def load_xml(self):

        ftmp = QFileDialog.getOpenFileName(self, 'Open file', \
            '%s/xml' %dpsd_dir, "xml files (*.xml)")
        if qt5:
            fxml = ftmp[0]
        else:
            fxml = str(ftmp)
        setup_d = xml.xml2dict(fxml)
        self.set_gui(setup_d)


    def save_xml(self):

        dpsd_dic = self.get_gui()
        ftmp = QFileDialog.getSaveFileName(self, 'Save file', \
            '%s/xml' %dpsd_dir, "xml files (*.xml)")
        if qt5:
            fxml = ftmp[0]
        else:
            fxml = str(ftmp)
        xml.dict2xml(dpsd_dic, fxml)


    def run(self):

        dpsd_dic = self.get_gui()
        self.dp = dpsd.DPSD(dpsd_dic)
        logger.info('Done calculation')


    def plot(self):

        if not hasattr(self, 'dp'):
            logger.error('Run code before plotting')
            return

        self.wid = plot_dpsd.plotWindow()
        fig1 = plot_dpsd.fig_pha(self.dp)
        fig2 = plot_dpsd.fig_phs(self.dp)
        fig3 = plot_dpsd.fig_cnt(self.dp)
        fig4 = plot_dpsd.fig_pmg(self.dp)
        fig5 = plot_dpsd.fig_win(self.dp)

        self.wid.addPlot('PH-PS separation'     , fig1)
        self.wid.addPlot('Pulse Height spectrum', fig2)
        self.wid.addPlot('Count rates'          , fig3)
        self.wid.addPlot('PM gain'              , fig4)
        self.wid.addPlot('Window length distribution', fig5)

        self.wid.show()
 

    def write_sf(self):

        if hasattr(self, 'dp'):
            self.dp.sfwrite(exp='git')
        else:
            logger.error('Run DPSD first, then write Shotfile')


if __name__ == '__main__':


    app = QApplication(sys.argv)
    main = DPSD()
    app.exec_()
