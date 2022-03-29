import os, sys, logging
import numba as nb

try:
    import Tkinter as tk
    import ttk
except:
    import tkinter as tk
    from tkinter import ttk
import dixm, dpsd
try:
    import aug_sfutils as sf
except:
    pass

fmt = logging.Formatter('%(asctime)s | %(name)s | %(levelname)s: %(message)s', '%H:%M:%S')
logger = logging.getLogger('DPSD_GUI')
hnd = logging.StreamHandler()
hnd.setFormatter(fmt)
logger.addHandler(hnd)
logger.setLevel(logging.INFO)

dpsd_dir = os.path.dirname(os.path.realpath(__file__))


class DPSD_GUI:


    def __init__(self):


        if __name__ == '__main__':
            topframe = tk.Tk()
            import ttk_style
        else:
            topframe = tk.Toplevel()

        topframe.title('DPSD')
        setup_en_d = dixm.DIX().xml2dict('/afs/ipp/home/g/git/DPSD/xml/default.xml')

        toolframe = ttk.Frame(topframe)
        entframe  = ttk.Frame(topframe)
        cbframe   = ttk.Frame(topframe)

        for frame in toolframe, entframe, cbframe:
            frame.pack(side=tk.TOP, fill =tk.BOTH, expand=1)

# Toolbar

        runfig  = tk.PhotoImage(file='%s/play.gif' %dpsd_dir)
        exitfig = tk.PhotoImage(file='%s/exit.gif' %dpsd_dir)
        plotfig = tk.PhotoImage(file='%s/plot.gif' %dpsd_dir)
        savefig = tk.PhotoImage(file='%s/save.gif' %dpsd_dir)
        btrun  = ttk.Button(toolframe, command=self.run, image=runfig)
        btexit = ttk.Button(toolframe, command=sys.exit, image=exitfig)
        btplot = ttk.Button(toolframe, command=self.plot, image=plotfig)
        btsave = ttk.Button(toolframe, command=self.write_sf, image=savefig)
        for but in btrun, btexit, btplot:
            but.pack(side=tk.LEFT)
        if 'aug_sfutils' in sys.modules:
            btsave.pack(side=tk.LEFT)

# Entries

        self.dpsd_d = {}

        cbdic = {'GraphPulse':False, 'GraphTime':False, 'GraphWin':False, 'SubtBaseline':True, 'LEDcorrection':True}

        jrow = 0
        entframe1 = ttk.Frame(entframe)
        entframe1.pack(side=tk.LEFT, padx='10 2')
        for key, val in setup_en_d.items():
            if key in cbdic.keys() or key == 'Slice':
                continue
            if jrow == 17:
                entframe1 = ttk.Frame(entframe)
                entframe1.pack(side=tk.LEFT, padx='10 2', anchor=tk.N)
                jrow = 0
            locframe = ttk.Frame(entframe1)
            locframe.pack(side=tk.TOP, anchor=tk.W)
            ttk.Label(locframe, text=key, width=12).pack(side=tk.LEFT, anchor=tk.W)
            if key == 'Path':
                self.dpsd_d[key] = ttk.Entry(locframe, width=32)
            else:
                self.dpsd_d[key] = ttk.Entry(locframe, width=6)
            self.dpsd_d[key].insert(0, val)
            self.dpsd_d[key].pack(side=tk.LEFT, anchor=tk.W)
            jrow += 1

# Checkbutton

        for key in ('SubtBaseline', 'LEDcorrection'):
            self.dpsd_d[key] = tk.BooleanVar()
            self.dpsd_d[key].set(cbdic[key])
            ttk.Checkbutton(cbframe, variable=self.dpsd_d[key], text=key).pack(side=tk.LEFT, padx='10 2', pady='10 2')

        topframe.mainloop()


    def run(self):

        dpsd_dic = {}
        for key, val in self.dpsd_d.items():
            dpsd_dic[key] = val.get()
        self.dp = dpsd.DPSD(dpsd_dic)


    def plot(self):

        if hasattr(self, 'dp'):
            self.dp.plot()
        else:
            logger.error('Run DPSD first, then plot')


    def write_sf(self):

        if hasattr(self, 'dp'):
            self.dp.sfwrite(exp='git')
        else:
            logger.error('Run DPSD first, then write Shotfile')


if __name__ == "__main__":

    DPSD_GUI()
