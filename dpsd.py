import sys, os, logging

import numpy as np
import matplotlib.pylab as plt
from matplotlib.patches import Rectangle
import aug_sfutils as sf
from aug_sfutils import sfhmod
import dixm, read_ha
from multiprocessing import Pool, cpu_count

ww = sf.WW()

fmt = logging.Formatter('%(asctime)s | %(name)s | %(levelname)s: %(message)s', '%H:%M:%S')
logger = logging.getLogger('DPSD')
hnd = logging.StreamHandler()
hnd.setFormatter(fmt)
logger.addHandler(hnd)
logger.setLevel(logging.INFO)

sig1d = ['neut1', 'neut2', 'gamma1', 'gamma2', 'led', 'pileup']


def BaselineCond(tuple_in):

    pulse, pulse_len, maxpos, bl_start, max_diff, aver1, max_sg, max_lg = tuple_in
    pulse_basestart = pulse_len - bl_start
    blstart_h = bl_start//2
    newpulse_len = 0
    if pulse_basestart >= maxpos:
        newpulse_len = pulse_len - blstart_h
    else:
        for j in range(maxpos, pulse_basestart):
            aver2 = np.average(pulse[j: j + bl_start])
            if np.abs(aver2 - aver1) < max_diff:
                if max_lg > j + blstart_h:
                    newpulse_len = max_lg
                else:
                    newpulse_len = j + blstart_h
                break
            if (j == pulse_basestart - 1) :
                newpulse_len = pulse_len - blstart_h
    short_int = np.trapz(pulse[maxpos: max_sg])
    long_int  = np.trapz(pulse[maxpos: max_lg])
    total_int = np.trapz(pulse[:newpulse_len])

    return np.array([short_int, long_int, total_int])


def PeakAlgo(nfront, ntail, nthres, pulse):

    pulse_width = nfront + ntail
    pulse_len = len(pulse)
    pile_out = 0
    jt = 0
    pulse_front = pulse[nfront : -ntail] - nthres
    pulse_max = np.maximum(pulse[:-pulse_width], pulse[pulse_width:])
    while jt < pulse_len-pulse_width:
        if pulse_front[jt] > pulse_max[jt]:
            pile_out += 1
            jt += pulse_width
        jt += 1
    return pile_out


class DPSD:


    def __init__(self, dic_in, sfw=True, exp='AUGD', fsfh='NSP00000.sfh'):

        
        self.nshot = int(dic_in['Acq1'])
        log = logging.getLogger('DPSD')

# Input

#        filepath = '%s/%s/%s' %(dic_in['Path'], str(dic_in['Acq1']), str(dic_in['Acq1']))
        shot100 = self.nshot//100
        filepath = '%s/%d/%d' %(dic_in['Path'], shot100, self.nshot)
        HAfile = '%s/HA_%d.dat' %(filepath, self.nshot)

        self.d_int = {}
        self.d_flt = {}
        int_list = ('BaselineStart', 'BaselineStart', 'BaselineEnd', \
            'ShortGate', 'LongGate', 'TotalGate', 'Threshold', 'Front', 'Tail', \
            'Marker', 'SaturationHigh', 'SaturationLow', 'ToFWindowLength', \
            'xChannels', 'yChannels', 'LineChange', \
            'LEDxmin', 'LEDxmax', 'LEDymin', 'LEDymax', 'LEDFront', 'LEDTail', 'LEDreference')
        flt_list = ('LEDdt', 'TimeBin', 'TBeg', 'TEnd', 'MaxDifference', 'Slope1', 'Slope2', 'Offset', )

        for key in int_list:
            self.d_int[key] = int(dic_in[key])
        for key in flt_list:
            self.d_flt[key] = float(dic_in[key])

        nxCh = self.d_int['xChannels']
        nyCh = self.d_int['yChannels']

        min_winlen = max(self.d_int['BaselineStart'], self.d_int['BaselineEnd'])
        ha = read_ha.READ_HA(HAfile, min_winlen=min_winlen, max_winlen=self.d_int['ToFWindowLength'])

        (tind, ) = np.where((ha.t_events >= self.d_flt['TBeg']) & (ha.t_events <= self.d_flt['TEnd']))

        time = ha.t_events[tind]
        print('TBeg = %8.4f' %time[0]) 
        print('TEnd = %8.4f' %time[-1]) 
        n_pulses = len(time)
        n_timebins = int((time[-1] - time[0])/self.d_flt['TimeBin'])
        n_led = int((time[-1] - time[0])/self.d_flt['LEDdt'])
        print('n_timebins = %d' %n_timebins)
        self.time = time[0] + self.d_flt['TimeBin']*(0.5 + np.arange(n_timebins))
        self.time_led = time[0] + self.d_flt['LEDdt']*(0.5 + np.arange(n_led))

        self.Xratio        = np.zeros(n_pulses, dtype=np.float32)
        self.Yratio        = np.zeros(n_pulses, dtype=np.float32)

        winlen = ha.winlen[tind]
        pulses = ha.pulses[tind]
        hist_len = np.bincount(winlen)

# Initialise
        flg_sat   = np.zeros(n_pulses)
        self.flg_peaks = np.zeros(n_pulses)

        self.pmgain = np.zeros(n_led, dtype=np.float32)

        LEDsumm = 0
        LEDaver = 0
        LEDcoeff = 0
        jtmark = 0

# =========== 
# Time loop
# =========== 

        logger.info('Starting time loop')

        print('# pulses: %6d' % n_pulses)
        ind1 = np.arange(self.d_int['BaselineStart'])
        maxpos = np.argmax(pulses, axis=1)
        pulse_max = np.max(pulses, axis=1)

        max_LG = np.minimum(maxpos + self.d_int['LongGate' ], winlen)
        max_SG = np.minimum(maxpos + self.d_int['ShortGate'], winlen)
        print(maxpos.shape)
        print(max_LG.shape)
        tof_win_len = self.d_int['ToFWindowLength']
        sat_high = float(dic_in['SaturationHigh'])
        sat_low  = float(dic_in['SaturationLow'])
        pulse_len = np.minimum(winlen, tof_win_len)
        pulse_baseend   = pulse_len - self.d_int['BaselineEnd']
        
        logger.info('Baseline subtraction')
        (ind_var, ) = np.where(pulse_len != tof_win_len)
        ind2 = np.arange(tof_win_len - self.d_int['BaselineEnd'], tof_win_len)
        ind_avg = np.unique(np.append(ind1, ind2))
        baseline = np.average(pulses[:, ind_avg], axis=1)

        for jpul in ind_var:
# Baseline subtraction
            ind2 = np.arange(pulse_baseend[jpul], pulse_len[jpul])
            ind_bl = np.unique(np.append(ind1, ind2))
            baseline[jpul] = np.average(pulses[jpul, ind_bl])

        pulses_flt = pulses.astype(np.float32) - baseline[:, None]

# Saturation detection
        logger.info('Saturation detection')
        (ind_sat_high, ) = np.where(np.max(pulses_flt, axis=1) > sat_high)
        (ind_sat_low, )  = np.where(np.min(pulses_flt, axis=1) < sat_low) 
        flg_sat[ind_sat_high] = 1
        flg_sat[ind_sat_low]  = 2

        logger.info('Baseline conditioned 2') 
        aver1 = np.average(pulses_flt[:, ind1], axis=1)

        timeout_pool = 120

        pool = Pool(cpu_count())
        out = pool.map_async(BaselineCond, [(pulses_flt[jpul], pulse_len[jpul], maxpos[jpul], self.d_int['BaselineStart'], self.d_flt['MaxDifference'], aver1[jpul], max_SG[jpul], max_LG[jpul]) for jpul in range(n_pulses)]).get(timeout_pool)

        out = np.array(out)
        self.ShortIntegral = out[: , 0]
        self.LongIntegral  = out[: , 1]
        self.TotalIntegral = out[: , 2]

        mark = np.float32(nxCh)/np.float32(self.d_int['Marker'])
        ind3 = np.where(self.LongIntegral > 0)[0]

# LED evaluation

        self.Xratio[ind3] = mark*self.TotalIntegral[ind3]
        self.Yratio[ind3] = np.float32(nyCh)*self.ShortIntegral[ind3]/self.LongIntegral[ind3];
        flg = {}

        flg['led'] =  \
            (self.Xratio > float(self.d_int['LEDxmin'])) & \
            (self.Xratio < float(self.d_int['LEDxmax'])) & \
            (self.Yratio > float(self.d_int['LEDymin'])) & \
            (self.Yratio < float(self.d_int['LEDymax'])) 

# Pile-up detection, use integer pulse (baseline just offsets whole pulse)

        logger.info('Pile-up detection')
        for jpul in range(n_pulses):

            if flg['led'][jpul]: #LED
                self.flg_peaks[jpul] = PeakAlgo( \
                    self.d_int['LEDFront'], self.d_int['LEDTail'], self.d_int['Threshold'], pulses[jpul])
            else:
                self.flg_peaks[jpul] = PeakAlgo( \
                    self.d_int['Front'], self.d_int['Tail'], self.d_int['Threshold'], pulses[jpul])

# LED correction

        logger.info('LED correction')
        self.cnt = {}
        led, tedges = np.histogram(time[flg['led']], bins=n_timebins, range=[self.time[0]-0.5*self.d_flt['TimeBin'], self.time[-1]+0.5*self.d_flt['TimeBin']])
        self.cnt['led'] = led.astype(np.float32)

        jtime_old = 0
        jtmark = 0
        LEDsumm = 0
        LEDamount = 0
        LEDcoeff = 0
        tled = ((time - time[0])/self.d_flt['LEDdt']).astype(np.int32) - 1
        for jpul in range(n_pulses):
            jtime = tled[jpul]
            if (jtime > jtime_old):
                if LEDamount > 0:
                    self.pmgain[jtime] = mark*np.float32(LEDsumm)/np.float32(LEDamount)
                    if LEDsumm > 0:
                        LEDcoeff = np.float32(self.d_int['LEDreference'])/self.pmgain[jtime]
                self.TotalIntegral[jtmark: jpul] *= LEDcoeff
                self.Xratio[jtmark: jpul] *= LEDcoeff
                jtmark = jpul
                LEDsumm = 0
                LEDamount = 0
                LEDcoeff = 0
            if flg['led'][jpul]: # LED single pulse
                LEDsumm += self.TotalIntegral[jpul]
                LEDamount += 1
            jtime_old = jtime

        flg_slope1 = (self.Xratio <= self.d_int['LineChange'])
        flg1n = (self.Yratio <= self.d_flt['Offset'] + self.d_flt['Slope1']*self.Xratio)
        offset2 = self.d_flt['Offset'] + self.d_flt['Slope1']*self.d_int['LineChange']
        flg2n = (self.Yratio <= offset2 + self.d_flt['Slope2']*(self.Xratio-self.d_int['LineChange']))

        flg1  =   flg_slope1  & flg1n
        flg2  = (~flg_slope1) & flg2n
        flg1g =   flg_slope1  & (~flg1n)
        flg2g = (~flg_slope1) & (~flg2n)
        flg['sat'] = (flg_sat > 0)
        flg['pileup'] = (self.flg_peaks > 1)
        flg['phys'] =  (~flg['sat']) & (~flg['led']) & (~flg['pileup'])
        flg['neut1']  = (flg1  + flg2 ) & (flg['phys'])
        flg['gamma1'] = (flg1g + flg2g) & (flg['phys'])

        cnt_list = ('neut1', 'gamma1', 'pileup', 'sat', 'phys')

        for spec in cnt_list:
            cnt, tedges = np.histogram(time[flg[spec]], bins=n_timebins, range=[self.time[0]-0.5*self.d_flt['TimeBin'], self.time[-1]+0.5*self.d_flt['TimeBin']])
            self.cnt[spec] = cnt.astype(np.float32)
            logger.info('%s %d', spec, np.sum(flg[spec]))
        logger.info('%s %d', 'LED', np.sum(flg['led']))

# Move to 1/ units
        for spec in self.cnt.keys():
            self.cnt[spec] /= self.d_flt['TimeBin']

        total = self.cnt['neut1'] + self.cnt['gamma1'] + self.cnt['led']
# Assuming pile-ups are all 2 events per window
        pup_frac = 1 + 2.*self.cnt['pileup']/total
        self.cnt['neut2' ] = pup_frac*self.cnt['neut1'] 
        self.cnt['gamma2'] = pup_frac*self.cnt['gamma1']


    def sfwrite(self, fsfh='NSP00000.sfh', exp='AUGD'):

        diag = 'NSP'
        nsp = sf.SFREAD('NSP', self.nshot, exp=exp)
        if nsp.status:
            logger.error('NSP shotfile for #%d exists already' %self.nshot)
            return

        nt = len(self.time)
        sfh = sfhmod.SFHMOD(fin=fsfh)
        for lbl in ['time'] + sig1d:
            sfh.modtime(lbl, nt)
        sfh.write(fout=fsfh)

        if ww.Open(exp, diag, self.nshot):
            status = ww.SetSignal('time', np.array(self.time, dtype=np.float32))
            for lbl in sig1d:
                print('Writing signal %s' %lbl)
                status = ww.SetSignal(lbl, np.array(self.cnt[lbl], dtype=np.float32))
            ww.Close()


    def plot(self):

# PHA run

        nbins = [4096, 1024]
        ranges = [[0, 4095], [0, 1023]]
        hpha, xedges, yedges = np.histogram2d(self.Xratio, self.Yratio, bins=nbins, range=ranges)
        hpha = np.rot90(hpha)
        hpha = np.flipud(hpha)
        Hmasked = np.ma.masked_where(hpha == 0, hpha) # Mask pixels with a value of zer#

        plt.figure('Pulse Shape separation', figsize = (13, 6.5))
        plt.pcolormesh(xedges, yedges, np.log10(Hmasked))
        cbar = plt.colorbar()

        xknot = self.d_int['LineChange']
        xline1 = [0, xknot]
        yknot = self.d_flt['Offset'] + self.d_flt['Slope1'] * (self.d_int['LineChange'] + 1)
        yline1 = [self.d_flt['Offset'], yknot]
        xline2 = [xknot, nbins[0]]
        yline2 = [yknot, yknot + self.d_flt['Slope2'] * (nbins[0] - xknot + 1)]
        plt.plot(xline1, yline1, 'r-')
        plt.plot(xline2, yline2, 'r-')
        cbar.ax.set_ylabel('Counts')

        xy = [self.d_int['LEDxmin'], self.d_int['LEDymin']]
        width  = self.d_int['LEDxmax'] - self.d_int['LEDxmin']
        height = self.d_int['LEDymax'] - self.d_int['LEDymin']
        led_box = Rectangle(xy, width, height, color='b', fill=False)
        plt.gca().add_patch(led_box)


        plt.figure('Count rates')
        for spec in sig1d:
            plt.plot(self.time, self.cnt[spec], label=spec)

        plt.legend()

        plt.figure('PM gain')
        plt.plot(self.time_led, self.pmgain/float(self.d_int['LEDreference']), 'r-')
        plt.show()


if __name__ == "__main__":

    if len(sys.argv) > 1:
        nshot = int(sys.argv[1])
    else:
        try:
            nshot = int(raw_input('Enter a shot number:'))
        except:
            print('Not a number')

    setup_en_d = dixm.DIX().xml2dict('/afs/ipp/home/g/git/DPSD/xml/default.xml')
    setup_en_d['TBeg'] = 0.
    setup_en_d['LEDdt'] = 0.01
    setup_en_d['TimeBin'] = 0.003
    setup_en_d['TEnd'] = 80.

#    dps = DPSD(setup_en_d, nshot, exp='git')
    dps = DPSD(setup_en_d, nshot, exp='AUGD')

# plots

    dps.plot()
