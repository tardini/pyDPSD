#!/usr/bin/env python

__author__  = 'Giovanni Tardini (Tel. +49 89 3299-1898)'
__version__ = '0.0.1'
__date__    = '29.03.2022'

import os, sys, logging, webbrowser

try:
    from PyQt5.QtWidgets import QMainWindow, QWidget, QApplication, QGridLayout, QMenu, QAction, QLabel, QPushButton, QLineEdit, QCheckBox, QFileDialog, QRadioButton, QButtonGroup, QTabWidget, QVBoxLayout
    from PyQt5.QtGui import QPixmap, QIcon, QIntValidator, QDoubleValidator
    from PyQt5.QtCore import Qt, QRect, QSize
    qt5 = True
except:
    from PyQt4.QtCore import Qt, QRect, QSize, QIntValidator, QDoubleValidator
    from PyQt4.QtGui import QPixmap, QIcon, QMainWindow, QWidget, QApplication, QGridLayout, QMenu, QAction, QLabel, QPushButton, QLineEdit, QCheckBox, QFileDialog, QRadioButton, QButtonGroup, QTabWidget, QVBoxLayout
    qt5 = False

import numpy as np
import dicxml, dpsd_run, plot_dpsd, plot_pulses
try:
    import aug_sfutils as sf
except:
    pass

os.environ['BROWSER'] = '/usr/bin/firefox'

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
        xmlMenu  = QMenu('&Setup', self)
        helpMenu = QMenu('&Help', self)
        menubar.addMenu(fileMenu)
        menubar.addMenu(xmlMenu)
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

        loadAction  = QAction('&Load Setup', xmlMenu)
        saveAction  = QAction('&Save Setup', xmlMenu)
        xmlMenu.addAction(loadAction)
        xmlMenu.addAction(saveAction)
        loadAction.triggered.connect(self.load_xml)
        saveAction.triggered.connect(self.save_xml)

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

        self.xml_d = dicxml.xml2dict('%s/xml/default.xml' %dpsd_dir)['main']
        self.setup_init = dicxml.xml2val_dic(self.xml_d)
        self.gui = {}
        for node in self.setup_init.keys():
            self.gui[node] = {}
        user = os.getenv('USER')

# Entry widgets

        self.rblists = {'SFexp':['AUGD', os.getenv('USER')]}

#----------
# I/O files
#----------

        jrow = 0
        node = 'io'

        key = 'HAfile'
        lbl = QLabel(self.xml_d[node][key]['@label'])
        self.gui[node][key] = QLineEdit(self.setup_init[node][key])
        input_layout.addWidget(lbl, jrow, 0)
        input_layout.addWidget(self.gui[node][key], jrow, 1, 1, 3)
        jrow += 1

        key = 'Shots'
        lbl = QLabel(self.xml_d[node][key]['@label'])
        self.gui[node][key] = QLineEdit(self.setup_init[node][key])
        input_layout.addWidget(lbl, jrow, 0)
        input_layout.addWidget(self.gui[node][key], jrow, 1, 1, 3)
        jrow += 1

        key = 'SFwrite'
        if 'aug_sfutils' in sys.modules:
            lbl = self.xml_d[node][key]['@label']
            self.gui[node][key] = QCheckBox(lbl)
            input_layout.addWidget(self.gui[node][key], jrow, 0, 1, 2)
            if self.setup_init[node][key]:
                self.gui[node][key].setChecked(True)
            jrow += 1

        key = 'SFforce'
        if 'aug_sfutils' in sys.modules:
            lbl = self.xml_d[node][key]['@label']
            self.gui[node][key] = QCheckBox(lbl)
            input_layout.addWidget(self.gui[node][key], jrow, 0, 1, 2)
            if self.setup_init[node][key]:
                self.gui[node][key].setChecked(True)
            jrow += 1

# Radiobutton

        key = 'SFexp'
        if 'aug_sfutils' in sys.modules:
            rblist = self.rblists[key]
            self.gui[node][key] = QButtonGroup(self)
            lbl = QLabel(self.xml_d[node][key]['@label'])
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

        entries = ['TimeBin', 'TBeg', 'TEnd', 'ToFWindowLength']
        self.new_tab(setup_layout, 'setup', entries=entries)

#-----
# Peak
#-----

        cb = ['SubtBaseline']
        entries = ['BaselineStart', 'BaselineEnd', 'Threshold', \
            'Front', 'Tail', 'SaturationHigh', 'SaturationLow', \
            'LongGate', 'ShortGate', 'MaxDifference']
        self.new_tab(peak_layout, 'peak', entries, checkbuts=cb)

#-----------
# Separation
#-----------

        entries = ['Marker', 'PH_nChannels', 'PS_nChannels', \
            'DDlower', 'DDupper', 'DTlower', 'DTupper', \
            'Slope1', 'Slope2', 'Offset', 'LineChange']
        self.new_tab(sep_layout, 'separation', entries)

#---------------
# LED correction
#---------------

        cb = ['LEDcorrection']
        entries = ['LEDdt', 'LEDFront', 'LEDTail', 'LEDreference', \
            'LEDxmin', 'LEDxmax', 'LEDymin', 'LEDymax']
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
            lbl = self.xml_d[node][key]['@label']
            self.gui[node][key] = QCheckBox(lbl)
            layout.addWidget(self.gui[node][key], jrow, 0, 1, 2)
            if self.setup_init[node][key]:
                self.gui[node][key].setChecked(True)
            jrow += 1

        for key in entries:
            val = self.setup_init[node][key]
            lbl = self.xml_d[node][key]['@label']
            qlbl = QLabel(lbl)
            qlbl.setFixedWidth(200)
            self.gui[node][key] = QLineEdit(str(val))
            self.gui[node][key].setFixedWidth(90)
            if isinstance(val, int):
                valid = QIntValidator()
                self.gui[node][key].setValidator(valid)
            elif isinstance(val, float):
                valid = QDoubleValidator()
                self.gui[node][key].setValidator(valid)
            layout.addWidget(qlbl         , jrow, 0)
            layout.addWidget(self.gui[node][key], jrow, 1)
            jrow += 1

        layout.setRowStretch(layout.rowCount(), 1)
        layout.setColumnStretch(layout.columnCount(), 1)


    def gui2xmld(self):
        '''Returns a dic of the xml-type, #text and @attr'''
        dpsd_dic = {}
        for node in self.gui.keys():
            dpsd_dic[node] = self.get_gui_tab(node)
        return dpsd_dic


    def get_gui_tab(self, node):

        node_dic = {}
        for key, val in self.gui[node].items():
            node_dic[key] = {}
            if isinstance(val, QLineEdit):
                node_dic[key]['#text'] = val.text()
            elif isinstance(val, QCheckBox):
                if val.isChecked():
                    node_dic[key]['#text'] = 'true'
                else:
                    node_dic[key]['#text'] = 'false'
            elif isinstance(val, QButtonGroup):
                bid = val.checkedId()
                node_dic[key]['#text'] = self.rblists[key][bid]
            node_dic[key]['@type' ] = self.xml_d[node][key]['@type']
            node_dic[key]['@label'] = self.xml_d[node][key]['@label']

        return node_dic


    def set_gui(self, xml_d):

        for node, val1 in xml_d.items():
            for key, vald in val1.items():
                if key not in self.gui[node].keys():
                    continue
                if '#text' in vald.keys():
                    val = vald['#text']
                else:
                    val = ''
                val = val.strip()
                val_low = val.lower()
                widget = self.gui[node][key]
                if isinstance(widget, QCheckBox):
                    if val_low == 'false':
                        widget.setChecked(False)
                    elif val_low == 'true':
                        widget.setChecked(True)
                elif isinstance(widget, QButtonGroup):
                    for but in widget.buttons():
                        if but.text().lower() == val_low:
                            but.setChecked(True)
                elif isinstance(widget, QLineEdit):
                    if val_low == '':
                        widget.setText(' ')
                    else:
                        widget.setText(val)
                elif isinstance(widget, QComboBox):
                    for index in range(widget.count()):
                        if widget.itemText(index).strip() == val.strip():
                            widget.setCurrentIndex(index)
                            break


    def load_xml(self):

        ftmp = QFileDialog.getOpenFileName(self, 'Open file', \
            '%s/xml' %dpsd_dir, "xml files (*.xml)")
        if qt5:
            fxml = ftmp[0]
        else:
            fxml = str(ftmp)
        setup_d = dicxml.xml2dict(fxml)
        self.set_gui(setup_d['main'])


    def save_xml(self):

        out_dic = {}
        out_dic['main'] = self.gui2xmld()
        ftmp = QFileDialog.getSaveFileName(self, 'Save file', \
            '%s/xml' %dpsd_dir, "xml files (*.xml)")
        if qt5:
            fxml = ftmp[0]
        else:
            fxml = str(ftmp)
        dicxml.dict2xml(out_dic, fxml)


    def run(self):

        dpsd_dic = dicxml.xml2val_dic(self.gui2xmld())
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
