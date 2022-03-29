import sys
import matplotlib.pylab as plt
import matplotlib
import numpy as np

from PyQt5.QtWidgets import QWidget, QTabWidget, QVBoxLayout
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar
from PyQt5.QtCore import QRect
matplotlib.rcParams['backend'] = 'Qt5Agg'
from matplotlib.patches import Rectangle

fsize   = 8
titsize = 10
lblsize = 10

sig1d = ['neut1', 'neut2', 'gamma1', 'gamma2', 'led', 'pileup']

class plotWindow(QWidget):


    def __init__(self):

        super().__init__()
        self.setGeometry(QRect(10, 10, 1200, 710))
        self.tabs = QTabWidget(self)

        self.setWindowTitle('Time traces')


    def addPlot(self, title, figure):

        new_tab = QWidget()
        layout = QVBoxLayout()
        new_tab.setLayout(layout)

        new_canvas = FigureCanvas(figure)
        new_toolbar = NavigationToolbar(new_canvas, new_tab)
#        new_canvas.mpl_connect('button_press_event', fconf.on_click)
        layout.addWidget(new_canvas)
        layout.addWidget(new_toolbar)
        self.tabs.addTab(new_tab, title)


def fig_pha(dpsd, color='#c00000'):

    fig = plt.figure(figsize=(12., 6.3), dpi=100)
    fig.subplots_adjust(left=0.05, bottom=0.07, right=0.98, top=0.94, hspace=0, wspace=0.4)
    fig.text(.5, .95, '#%d' %dpsd.nshot, ha='center')

    nbins = [dpsd.d_int['xChannels'], dpsd.d_int['yChannels']]
    ranges = [[-0.5, nbins[0]+0.5], [-0.5, nbins[1]+0.5]]
    hpha, xedges, yedges = np.histogram2d(dpsd.Xratio, dpsd.Yratio, bins=nbins, range=ranges)
    hpha = np.rot90(hpha)
    hpha = np.flipud(hpha)
    Hmasked = np.ma.masked_where(hpha == 0, hpha) # Mask pixels with a value of zer#

    plt.pcolormesh(xedges, yedges, np.log10(Hmasked))
    cbar = plt.colorbar()

    xknot = dpsd.d_int['LineChange']
    xline1 = [0, xknot]
    yknot = dpsd.d_flt['Offset'] + dpsd.d_flt['Slope1'] * (dpsd.d_int['LineChange'] + 1)
    yline1 = [dpsd.d_flt['Offset'], yknot]
    xline2 = [xknot, nbins[0]]
    yline2 = [yknot, yknot + dpsd.d_flt['Slope2'] * (nbins[0] - xknot + 1)]
    plt.plot(xline1, yline1, 'r-')
    plt.plot(xline2, yline2, 'r-')

    xy = [dpsd.d_int['LEDxmin'], dpsd.d_int['LEDymin']]
    width  = dpsd.d_int['LEDxmax'] - dpsd.d_int['LEDxmin']
    height = dpsd.d_int['LEDymax'] - dpsd.d_int['LEDymin']
    led_box = Rectangle(xy, width, height, color='b', fill=False)
    plt.gca().add_patch(led_box)

    return fig


def fig_phs(dpsd, color='#c00000', ymax=2, titles=None):

    fig = plt.figure(figsize=(8.0, 6.3), dpi=100)

    fig.subplots_adjust(left=0.05, bottom=0.08, right=0.98, top=0.92, hspace=0, wspace=0.28)
    fig.text(.5, .95, '#%d' %dpsd.nshot, ha='center')

    ymax = 0
    for spec in ['neut1', 'gamma1', 'led']:
        plt.plot(dpsd.phs[spec], label=spec)
        ymax = max(ymax, np.max(dpsd.phs[spec][1:]))
    plt.ylim([0, ymax])
    plt.legend()

    return fig


def fig_cnt(dpsd, color='#c00000', ymax=2, titles=None):

    fig = plt.figure(figsize=(8.0, 6.3), dpi=100)

    fig.subplots_adjust(left=0.05, bottom=0.08, right=0.98, top=0.92, hspace=0, wspace=0.28)
    fig.text(.5, .95, '#%d' %dpsd.nshot, ha='center')

    for spec in sig1d:
        plt.plot(dpsd.time, dpsd.cnt[spec], label=spec)
    plt.legend()

    return fig


def fig_pmg(dpsd):

    fig = plt.figure(figsize=(8.0, 6.3), dpi=100)

    fig.subplots_adjust(left=0.05, bottom=0.08, right=0.98, top=0.92, hspace=0, wspace=0.28)
    fig.text(.5, .95, '#%d' %dpsd.nshot, ha='center')

    plt.plot(dpsd.time_led, dpsd.pmgain/float(dpsd.d_int['LEDreference']), 'r-')

    return fig
