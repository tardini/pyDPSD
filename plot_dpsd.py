import sys
import matplotlib
import numpy as np

try:
    from PyQt5.QtWidgets import QWidget, QTabWidget, QVBoxLayout
    from PyQt5.QtCore import QRect
    matplotlib.use('Qt5Agg')
    from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
    from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar
except:
    from PyQt4.QtGui import QWidget, QTabWidget, QVBoxLayout
    from PyQt4.QtCore import QRect
    matplotlib.use('Qt4Agg')
    from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas
    from matplotlib.backends.backend_qt4agg import NavigationToolbar2QT as NavigationToolbar

import matplotlib.pylab as plt
from matplotlib.patches import Rectangle

fsize   = 8
titsize = 10
lblsize = 10
fig_size = (8.2, 6.05)


class plotWindow(QWidget):


    def __init__(self):

        if sys.version_info[0] == 3:
            super().__init__()
        else:
            super(QWidget, self).__init__()
        self.setGeometry(QRect(80, 30, 850, 705))
        self.tabs = QTabWidget(self)
        self.tabs.setStyleSheet("QTabBar::tab { width: 165 }")
        self.setWindowTitle('DPSD output')


    def addPlot(self, title, figure):

        new_tab = QWidget()
        layout = QVBoxLayout()
        new_tab.setLayout(layout)

        new_canvas = FigureCanvas(figure)
        new_toolbar = NavigationToolbar(new_canvas, new_tab)
        layout.addWidget(new_canvas)
        layout.addWidget(new_toolbar)
        self.tabs.addTab(new_tab, title)


def fig_pha(dpsd, color='#c00000'):

    sep_d = dpsd.setup['separation']
    led_d = dpsd.setup['led']
    fig_pha = plt.figure(figsize=fig_size, dpi=100)
    fig_pha.subplots_adjust(left=0.1, bottom=0.1, right=0.98, top=0.94, hspace=0, wspace=0.4)
    if hasattr(dpsd, 'nshot'):
        fig_pha.text(.5, .95, '#%d' %dpsd.nshot, ha='center')

    nbins = [sep_d['PH_nChannels'], sep_d['PS_nChannels']]
    ranges = [[-0.5, nbins[0]+0.5], [-0.5, nbins[1]+0.5]]
    hpha, xedges, yedges = np.histogram2d(dpsd.PulseHeight, dpsd.PulseShape, bins=nbins, range=ranges)
    hpha = np.flipud(np.rot90(hpha))
    hpha[hpha == 0] = np.nan
    plt.xlim([0, nbins[0]])
    plt.ylim([0, nbins[1]])
    plt.pcolormesh(xedges, yedges, np.log10(hpha), cmap=matplotlib.cm.jet)
    cbar = plt.colorbar()

    xknot = sep_d['LineChange']
    xline1 = [0, xknot]
    yknot = sep_d['Offset'] + sep_d['Slope1'] * (sep_d['LineChange'] + 1)
    yline1 = [sep_d['Offset'], yknot]
    xline2 = [xknot, nbins[0]]
    yline2 = [yknot, yknot + sep_d['Slope2'] * (nbins[0] - xknot + 1)]
    plt.plot(xline1, yline1, 'r-')
    plt.plot(xline2, yline2, 'r-')
    for lbl in ('DDlower', 'DDupper'):
        plt.plot([sep_d[lbl], sep_d[lbl]], [0, nbins[1]], 'g-')
    for lbl in ('DTlower', 'DTupper'):
        plt.plot([sep_d[lbl], sep_d[lbl]], [0, nbins[1]], 'm-')

    xy = [led_d['LEDxmin'], led_d['LEDymin']]
    width  = led_d['LEDxmax'] - led_d['LEDxmin']
    height = led_d['LEDymax'] - led_d['LEDymin']
    led_box = Rectangle(xy, width, height, color='b', fill=False)
    plt.gca().add_patch(led_box)
    plt.xlabel('Pulse Height')
    plt.ylabel('Pulse Shape')

    return fig_pha


def fig_phs(dpsd, color='#c00000', ymax=2, titles=None):

    fig_phs = plt.figure(figsize=fig_size, dpi=100)

    fig_phs.subplots_adjust(left=0.1, bottom=0.1, right=0.98, top=0.92, hspace=0, wspace=0.28)
    fig_phs.text(.5, .95, '#%d' %dpsd.nshot, ha='center')

    ymax = 0
    for spec in ['neut1', 'gamma1', 'led', 'DT']:
        plt.plot(dpsd.phs[spec], label=spec)
        ymax = max(ymax, np.max(dpsd.phs[spec][1:]))
    plt.xlim([0, dpsd.setup['separation']['PH_nChannels']])
    plt.ylim([0, ymax])
    plt.xlabel('Pulse Height')
    plt.ylabel('Occurrences')
    plt.legend()

    return fig_phs


def fig_cnt(dpsd, color='#c00000', ymax=2, titles=None):

    fig_cnt = plt.figure(figsize=fig_size, dpi=100)

    fig_cnt.subplots_adjust(left=0.1, bottom=0.1, right=0.98, top=0.92, hspace=0, wspace=0.28)
    if hasattr(dpsd, 'nshot'):
        fig_cnt.text(.5, .95, '#%d' %dpsd.nshot, ha='center')

    ymax = 0
    for spec in ['neut1', 'neut2', 'gamma1', 'gamma2', 'led', 'pileup', 'DD', 'DT']:
        plt.plot(dpsd.time_cnt, dpsd.cnt[spec], label=spec)
        ymax = max(ymax, np.max(dpsd.cnt[spec]))
    plt.xlim([dpsd.time_cnt[0], dpsd.time_cnt[-1]])
    plt.ylim([0, ymax])
    plt.xlabel('Time [s]')
    plt.ylabel('Count rate [1/s]')
    plt.legend()

    return fig_cnt


def fig_pmg(dpsd):

    fig_pmg = plt.figure(figsize=fig_size, dpi=100)

    fig_pmg.subplots_adjust(left=0.1, bottom=0.1, right=0.98, top=0.92, hspace=0, wspace=0.28)
    if hasattr(dpsd, 'nshot'):
        fig_pmg.text(.5, .95, '#%d' %dpsd.nshot, ha='center')

    plt.plot(dpsd.time_led, dpsd.pmgain/float(dpsd.setup['led']['LEDreference']), 'r-')
    plt.xlim([dpsd.time_led[0], dpsd.time_led[-1]])
    plt.ylim([0, 1.5])
    plt.xlabel('Time [s]')
    plt.ylabel('Photomultiplier gain')
    return fig_pmg


def fig_win(dpsd):

    fig_win = plt.figure(figsize=fig_size, dpi=100)

    fig_win.subplots_adjust(left=0.1, bottom=0.1, right=0.98, top=0.92, hspace=0, wspace=0.28)
    if hasattr(dpsd, 'nshot'):
        fig_win.text(.5, .95, '#%d' %dpsd.nshot, ha='center')

    win_min = np.min(dpsd.winlen)
    win_max = np.max(dpsd.winlen)
    plt.hist(dpsd.winlen, bins=win_max)
    plt.xlim([win_min, win_max])
    plt.xlabel('Window length [#samples]')
    plt.ylabel('Occurrences')
    return fig_win
