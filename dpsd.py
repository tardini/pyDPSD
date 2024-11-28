#!/usr/bin/env python

__author__  = 'Giovanni Tardini (Tel. +49 89 3299-1898)'
__version__ = '0.0.1'
__date__    = '29.03.2022'

import os, sys, logging, webbrowser, json

try:
    from PyQt5.QtWidgets import QMainWindow, QWidget, QApplication, QGridLayout, QMenu, QAction, QLabel, QPushButton, QLineEdit, QCheckBox, QSpinBox, QDoubleSpinBox, QFileDialog, QRadioButton, QButtonGroup, QTabWidget, QVBoxLayout
    from PyQt5.QtGui import QPixmap, QIcon
    from PyQt5.QtCore import Qt, QRect, QSize, QLocale
    qt5 = True
except:
    from PyQt4.QtCore import Qt, QRect, QSize, QLocale
    from PyQt4.QtGui import QPixmap, QIcon, QMainWindow, QWidget, QApplication, QGridLayout, QMenu, QAction, QLabel, QPushButton, QLineEdit, QCheckBox, QSpinBox, QDoubleSpinBox, QFileDialog, QRadioButton, QButtonGroup, QTabWidget, QVBoxLayout
    qt5 = False

import numpy as np
import dpsd_run, plot_dpsd, plot_pulses
try:
    import aug_sfutils as sf
except:
    pass

os.environ['BROWSER'] = '/usr/bin/firefox'

usLocale = QLocale('us')

fmt = logging.Formatter('%(asctime)s | %(name)s | %(levelname)s: %(message)s', '%H:%M:%S')
hnd = logging.StreamHandler()
hnd.setFormatter(fmt)
logger = logging.getLogger('DPSD_GUI')
logger.addHandler(hnd)
logger.setLevel(logging.INFO)

dpsd_dir = os.path.dirname(os.path.realpath(__file__))


class DPSD(QMainWindow):


    def __init__(self):

        if sys.version_info[0] == 3:
            super().__init__()
        else:
            super(QMainWindow, self).__init__()

        self.setLocale(usLocale)

        xwin  = 603
        yhead = 44
        yline = 30
        ybar  = 48
        ywin  = xwin + yhead + ybar

        qhead  = QWidget(self)
        qbar   = QWidget(self)
        qtabs  = QTabWidget(self)
        qhead.setGeometry(QRect(0,     0, xwin, yhead))
        qbar.setGeometry(QRect(0, yhead, xwin, ybar))
        qtabs.setGeometry(QRect(0, yhead+ybar, xwin, ywin-yhead-ybar))
        qtabs.setStyleSheet("QTabBar::tab { width: 120 }")
        header_grid = QGridLayout(qhead) 
        tbar_grid   = QGridLayout(qbar) 

#-----
# Tabs
#-----

        qinput = QWidget()
        input_layout = QGridLayout()
        qinput.setLayout(input_layout)
        qtabs.addTab(qinput, 'I/O files')

        qsetup = QWidget()
        setup_layout = QGridLayout()
        qsetup.setLayout(setup_layout)
        qtabs.addTab(qsetup, 'Setup')

        qpeak = QWidget()
        peak_layout = QGridLayout()
        qpeak.setLayout(peak_layout)
        qtabs.addTab(qpeak, 'Peak detection')

        qsep = QWidget()
        sep_layout = QGridLayout()
        qsep.setLayout(sep_layout)
        qtabs.addTab(qsep, 'n-g separation')

        qled = QWidget()
        led_layout = QGridLayout()
        qled.setLayout(led_layout)
        qtabs.addTab(qled, 'LED correction')

#--------
# Menubar
#--------

        menubar = self.menuBar()
        fileMenu = QMenu('&File', self)
        jsonMenu  = QMenu('&Setup', self)
        helpMenu = QMenu('&Help', self)
        menubar.addMenu(fileMenu)
        menubar.addMenu(jsonMenu)
        menubar.addMenu(helpMenu)

        runAction  = QAction('&Run'     , fileMenu)
        plotAction = QAction('&Plots'   , fileMenu)
        ppulAction = QAction('&Pulse analysis', fileMenu)
        wsfAction  = QAction('&Write SF', fileMenu)
        exitAction = QAction('&Exit'    , fileMenu)
        runAction.triggered.connect(self.run)
        plotAction.triggered.connect(self.plot)
        ppulAction.triggered.connect(self.plot_pulses)
        wsfAction.triggered.connect(self.write_sf)
        exitAction.triggered.connect(sys.exit)
        fileMenu.addAction(runAction)
        fileMenu.addAction(plotAction)
        fileMenu.addAction(ppulAction)
        if 'aug_sfutils' in sys.modules:
            fileMenu.addAction(wsfAction)
        fileMenu.addSeparator()
        fileMenu.addAction(exitAction)

        loadAction  = QAction('&Load Setup', jsonMenu)
        saveAction  = QAction('&Save Setup', jsonMenu)
        jsonMenu.addAction(loadAction)
        jsonMenu.addAction(saveAction)
        loadAction.triggered.connect(self.load_json)
        saveAction.triggered.connect(self.save_json)

        aboutAction = QAction('&Web docu', helpMenu)
        aboutAction.triggered.connect(self.about)
        helpMenu.addAction(aboutAction)

        header_grid.addWidget(menubar, 0, 0, 1, 10)

# Icons
        dum_lbl  = QLabel(200*' ')
        fmap = {'exec': self.run, 'plot': self.plot, 'pulse': self.plot_pulses, \
            'save': self.write_sf, 'exit': sys.exit}
        for jpos, lbl in enumerate(['exec', 'exit', 'plot', 'pulse']):
            but = QPushButton()
            but.setIcon(QIcon('%s/%s.gif' %(dpsd_dir, lbl)))
            but.setIconSize(QSize(ybar, ybar))
            but.clicked.connect(fmap[lbl])
            tbar_grid.addWidget(but, 0, jpos)
        jpos += 1
        if 'aug_sfutils' in sys.modules:
            lbl = 'save'
            but = QPushButton()
            but.setIcon(QIcon('%s/%s.gif' %(dpsd_dir, lbl)))
            but.setIconSize(QSize(ybar, ybar))
            but.clicked.connect(fmap[lbl])
            tbar_grid.addWidget(but, 0, jpos)
            jpos += 1
        tbar_grid.addWidget(dum_lbl,  0, jpos, 1, 10)

# User options

        f_json = '%s/settings/default.json' %dpsd_dir
        with open(f_json) as fjson:
            self.setup_init = json.load(fjson)

        self.gui = {}
        for node in self.setup_init.keys():
            self.gui[node] = {}
        user = os.getenv('USER')

# Entry widgets

        self.rblists = {'Shotfile exp':['AUGD', os.getenv('USER')]}

#----------
# I/O files
#----------

        jrow = 0
        node = 'io'

        for key in ('HA*.dat file', 'Shots'):
            lbl = QLabel(key)
            self.gui[node][key] = QLineEdit(self.setup_init[node][key])
            input_layout.addWidget(lbl, jrow, 0)
            input_layout.addWidget(self.gui[node][key], jrow, 1, 1, 3)
            jrow += 1

        key = 'Write shotfiles'
        if 'aug_sfutils' in sys.modules:
            self.gui[node][key] = QCheckBox(key)
            input_layout.addWidget(self.gui[node][key], jrow, 0, 1, 2)
            if self.setup_init[node][key]:
                self.gui[node][key].setChecked(True)
            jrow += 1

        key = 'Force SF write'
        if 'aug_sfutils' in sys.modules:
            self.gui[node][key] = QCheckBox(key)
            input_layout.addWidget(self.gui[node][key], jrow, 0, 1, 2)
            if self.setup_init[node][key]:
                self.gui[node][key].setChecked(True)
            jrow += 1

# Radiobutton

        key = 'Shotfile exp'
        if 'aug_sfutils' in sys.modules:
            rblist = self.rblists[key]
            self.gui[node][key] = QButtonGroup(self)
            lbl = QLabel(key)
            input_layout.addWidget(lbl, jrow, 0)
            for jcol, val in enumerate(rblist):
                but = QRadioButton(val)
                if jcol == 0:
                    but.setChecked(True)
                if self.setup_init[node][key].lower() == val.lower():
                    but.setChecked(True)
                input_layout.addWidget(but, jrow, jcol+1)
                self.gui[node][key].addButton(but)
                self.gui[node][key].setId(but, jcol)
            jrow += 1

        input_layout.setRowStretch(input_layout.rowCount(), 1)

#------
# Setup
#------

        entries = ['Time step', 'Start time', 'End time', '#samples for analysis']
        self.new_tab(setup_layout, 'setup', entries=entries)

#-----
# Peak
#-----

        cb = ['Subtract baseline']
        entries = ['Baseline start', 'Baseline end', 'Threshold', \
            'Front', 'Tail', 'Saturation upper limit', 'Saturation lower limit', \
            'Long gate', 'Short gate', 'Maximum difference']
        self.new_tab(peak_layout, 'peak', entries, checkbuts=cb)

#-----------
# Separation
#-----------

        entries = ['Marker', '#bins Pulse Height', '#bins Pulse Shape', \
            'Lower PH-limit for DD', 'Upper PH-limit for DD', \
            'Lower PH-limit for DT', 'Upper PH-limit for DT', \
            'Slope of 1st sep.line', 'Slope of 2nd sep.line', 'Offset of 1st sep.line', \
            'Bin line1 -> line2']
        self.new_tab(sep_layout, 'separation', entries)

#---------------
# LED correction
#---------------

        cb = ['LED correction']
        entries = ['LED time sampling', 'LED front', 'LED tail', 'LED reference bin', \
                   'Min PS bin for LED detection', 'Max PS bin for LED detection', \
                   'Min PH bin for LED detection', 'Max PH bin for LED detection']
        self.new_tab(led_layout, 'led', entries, checkbuts=cb)

        self.setStyleSheet("QLabel { width: 4 }")
        self.setStyleSheet("QLineEdit { width: 4 }")
        self.setGeometry(10, 10, xwin, ywin)
        self.setWindowTitle('DPSD')
        self.show()


    def about(self):

        webbrowser.open('https://www.aug.ipp.mpg.de/~git/ne213/dpsd/index.html')


    def new_tab(self, layout, node, entries=[], checkbuts=[]):
# Checkbutton

        jrow = 0
        for key in checkbuts:
            self.gui[node][key] = QCheckBox(key)
            layout.addWidget(self.gui[node][key], jrow, 0, 1, 2)
            if self.setup_init[node][key]:
                self.gui[node][key].setChecked(True)
            jrow += 1

        for key in entries:
            val = self.setup_init[node][key]
            qlbl = QLabel(key)
            qlbl.setFixedWidth(200)
            self.gui[node][key] = QLineEdit(str(val))
            self.gui[node][key].setFixedWidth(90)
            if isinstance(val, int):
                self.gui[node][key] = QSpinBox()
                self.gui[node][key].setRange(-int(1e3), int(1e7))
                self.gui[node][key].setValue(val)
            elif isinstance(val, float):
                self.gui[node][key] = QDoubleSpinBox()
                self.gui[node][key].setRange(-1000., 1000.)
                self.gui[node][key].setDecimals(4)
                self.gui[node][key].setValue(val)
            layout.addWidget(qlbl         , jrow, 0)
            layout.addWidget(self.gui[node][key], jrow, 1)
            jrow += 1

        layout.setRowStretch(layout.rowCount(), 1)
        layout.setColumnStretch(layout.columnCount(), 1)


    def gui2json(self):

        json_d = {}
        for node, val in self.gui.items():
            json_d[node] = self.get_gui_tab(node)
        return json_d


    def get_gui_tab(self, node):

        node_dic = {}
        for key, val in self.gui[node].items():
            node_dic[key] = {}
            if isinstance(val, QLineEdit):
                node_dic[key] = val.text()
            elif isinstance(val, QSpinBox):
                node_dic[key] = int(val.text())
            elif isinstance(val, QDoubleSpinBox):
                node_dic[key] = float(val.text())
            elif isinstance(val, QCheckBox):
                node_dic[key] = val.isChecked()
            elif isinstance(val, QButtonGroup):
                bid = val.checkedId()
                node_dic[key] = self.rblists[key][bid]

        return node_dic


    def set_gui(self, json_d):

        for node, val1 in json_d.items():
            for key, vald in val1.items():
                if key not in self.gui[node].keys():
                    continue
                widget = self.gui[node][key]
                if isinstance(widget, QCheckBox):
                    widget.setChecked(vald)
                elif isinstance(widget, QButtonGroup):
                    for but in widget.buttons():
                        if but.text() == vald:
                            but.setChecked(True)
                elif isinstance(widget, QLineEdit):
                    if vald:
                        widget.setText(str(vald))
                elif isinstance(widget, QComboBox):
                    for index in range(widget.count()):
                        if widget.itemText(index).strip() == vald.strip():
                            widget.setCurrentIndex(index)
                            break


    def load_json(self):

        ftmp = QFileDialog.getOpenFileName(self, 'Open file', \
            '%s/settings' %dpsd_dir, "json files (*.json)")
        if qt5:
            f_json = ftmp[0]
        else:
            f_json = str(ftmp)

        with open(f_json) as fjson:
            setup_d = json.load(fjson)
        self.set_gui(setup_d)


    def save_json(self):

        out_dic = self.gui2json()
        ftmp = QFileDialog.getSaveFileName(self, 'Save file', \
            '%s/settings' %dpsd_dir, "json files (*.json)")
        if qt5:
            f_json = ftmp[0]
        else:
            f_json = str(ftmp)
        with open(f_json, 'w') as fjson:
            fjson.write(json.dumps(out_dic))


    def run(self):

        dpsd_dic = self.gui2json()
        self.dp = dpsd_run.DPSD(dpsd_dic)
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
        self.wid.addPlot('Pulse Height spectra', fig2)
        self.wid.addPlot('Count rates'          , fig3)
        self.wid.addPlot('PM gain'              , fig4)
        self.wid.addPlot('Window lengths', fig5)

        self.wid.show()

 
    def plot_pulses(self):

        if not hasattr(self, 'dp'):
            logger.error('Run code before plotting')
            return

        self.pul = plot_pulses.plotWindow(self.dp)
        self.pul.show()
 

    def write_sf(self):

        if hasattr(self, 'dp'):
            dic = self.get_gui_tab('io')
            if self.dp.status:
                self.dp.sfwrite(exp=dic['SFexp'], force=dic['SFforce'])
        else:
            logger.error('Run DPSD first, then write Shotfile')


if __name__ == '__main__':


    app = QApplication(sys.argv)
    main = DPSD()
    app.exec_()
