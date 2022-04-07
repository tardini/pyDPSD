import os, sys
import matplotlib
import numpy as np

try:
    from PyQt5.QtWidgets import QWidget, QTabWidget, QVBoxLayout, QGridLayout, QPushButton, QLabel, QProgressBar, QSlider
    from PyQt5.QtGui import QPixmap, QIcon
    from PyQt5.QtCore import Qt, QRect, QSize
    matplotlib.use('Qt5Agg')
    from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
    from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar
except:
    from PyQt4.QtGui import QPixmap, QIcon, QWidget, QTabWidget, QVBoxLayout, QGridLayout, QPushButton, QLabel, QProgressBar, QSlider
    from PyQt4.QtCore import Qt, QRect, QSize
    matplotlib.use('Qt4Agg')
    from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas
    from matplotlib.backends.backend_qt4agg import NavigationToolbar2QT as NavigationToolbar

import matplotlib.pylab as plt
from matplotlib.patches import Rectangle
from matplotlib.figure import Figure

fsize   = 8
titsize = 10
lblsize = 10

dpsd_dir = os.path.dirname(os.path.realpath(__file__))



class plotWindow(QWidget):


    def __init__(self, dpsd):

        self.dp = dpsd

        print(matplotlib.rcParams['backend'])
        if matplotlib.rcParams['backend'] == 'Qt5Agg':
            super().__init__()
        else:
            super(QWidget, self).__init__()
        xwin = 900
        ywin = 710
        self.timeout = 1.e-10
        xicon = 40
        yicon = 50
        ybar = 20

        self.setGeometry(QRect(80, 30, xwin, ywin))

        self.setWindowTitle('DPSD pulse analysis')

        layout = QVBoxLayout()
        self.setLayout(layout)

        self.fig_pul = Figure(figsize=(5., 4.3), dpi=100)
        self.canvas = FigureCanvas(self.fig_pul)
        new_toolbar = NavigationToolbar(self.canvas, self)
        tbar = QWidget(self)
        tbar_grid = QGridLayout(tbar)
        slider_layout = QVBoxLayout()
        self.progress = QLabel('0%')
        self.progress.setAlignment(Qt.AlignCenter)
        self.progress.setFixedHeight(ybar)
        self.slider = QSlider(Qt.Horizontal)
        self.slider.sliderReleased.connect(self.get_slvalue)
        slider_layout.addWidget(self.progress)
        slider_layout.addWidget(self.slider)

        tbar.setFixedHeight(yicon)

# Icons
        dum_lbl  = QLabel(200*' ')
        fmap = {'play': self.play, 'forward': self.forward, \
            'backward': self.backward, 'pause': self.pause}

        self.but = {}
        for jpos, lbl in enumerate(['backward', 'play', 'forward']):
            self.but[lbl] = QPushButton()
            self.but[lbl].setIcon(QIcon('%s/%s.gif' %(dpsd_dir, lbl)))
            self.but[lbl].setIconSize(QSize(xicon, yicon))
            self.but[lbl].clicked.connect(fmap[lbl])
            tbar_grid.addWidget(self.but[lbl], 0, jpos)
        tbar_grid.addWidget(dum_lbl, 0, jpos+1)

        layout.addWidget(self.canvas)
        layout.addLayout(slider_layout)
        layout.addWidget(tbar)
        layout.addWidget(new_toolbar)

#        self.canvas.mpl_connect('button_press_event', self.pause)
        self.ax = {}
        self.line = {}
        self.nt, self.xlen = self.dp.pulses.shape
        self.slider.setMaximum(self.nt)
        for j in range(4):
            self.ax[j] = self.fig_pul.add_subplot(2, 2, j + 1)
            self.ax[j].set_xlim([0, self.xlen])
            self.line[j], = self.ax[j].plot(np.arange(self.xlen), np.zeros(self.xlen))
        self.ax[0].set_title('Neutrons')
        self.ax[1].set_title('Gamma')
        self.ax[2].set_title('Pile-up')
        self.ax[3].set_title('LED')

        self.ftext = self.fig_pul.text(0.5, 0.95, 'Time: %7.5f' %self.dp.time[0], ha='center')

        self.jt = 0
        self.stop = True
        self.update_plot()


    def pause(self):
        self.but['play'].setIcon(QIcon('%s/play.gif' %dpsd_dir))
        self.but['play'].disconnect()
        self.but['play'].clicked.connect(self.play)
        self.stop = True

    def play(self):
        self.stop = False
        self.but['play'].setIcon(QIcon('%s/pause.gif' %dpsd_dir))
        self.but['play'].disconnect()
        self.but['play'].clicked.connect(self.pause)
        while self.jt < self.nt-1 and not self.stop:
            self.jt += 1
            self.update_plot()
            self.canvas.start_event_loop(self.timeout)

    def forward(self):
        self.jt += 1
        self.jt = self.jt%self.nt
        self.update_plot()
        self.pause()

    def backward(self):
        self.jt -= 1
        self.jt = self.jt%self.nt
        self.update_plot()
        self.pause()

    def update_plot(self):
        jplot = self.dp.event_type[self.jt]
        if jplot < 0:
            return
        self.progress.setText('%d%%' %((self.jt*100)//self.nt))
        self.slider.setValue(self.jt)
        pulse = self.dp.pulses[self.jt]
        self.line[jplot].set_ydata(pulse)
        self.ax[jplot].set_ylim([0, np.max(pulse)])
        self.ftext.set_text('Time=%7.5f' %self.dp.time[self.jt])
        self.canvas.draw()

    def get_slvalue(self):
        self.jt = self.slider.value()
        self.update_plot()
        self.pause()
