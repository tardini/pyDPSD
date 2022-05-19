import sys, os, logging

import numpy as np
import numba as nb
import dixm, read_ha


fmt = logging.Formatter('%(asctime)s | %(name)s | %(levelname)s: %(message)s', '%H:%M:%S')
logger = logging.getLogger('DPSD')
hnd = logging.StreamHandler()
hnd.setFormatter(fmt)
logger.addHandler(hnd)
logger.setLevel(logging.DEBUG)

sig1d = ['neut1', 'neut2', 'gamma1', 'gamma2', 'led', 'pileup']
dpsd_dir = os.path.dirname(os.path.realpath(__file__))

@nb.njit
def slice_trapz(a, bnd_l, bnd_r):
    b = np.zeros(a.shape[0])
    for j in range(a.shape[0]):
        for i in range(bnd_l[j]+1, bnd_r[j]-1):
            b[j] += a[j, i]
        b[j] += 0.5*(a[j, bnd_l[j]] + a[j, bnd_r[j]-1])
    return b


@nb.njit
def led_correction(dtled, dxCh, led_ref, time, totalintegral, flg_led):

    jled_old = 0
    jtmark = 0
    LEDsumm = 0
    LEDamount = 0
    LEDcoeff = 0
    tled = ((time - time[0])/dtled).astype(np.int32)
    n_pulses = len(time)
    n_led = int((time[-1] - time[0])/dtled)
    pulseheight = dxCh*totalintegral
    pmgain = np.zeros(n_led, dtype=np.float32)

    for jpul in range(n_pulses):
        jled = tled[jpul]
        if jled > jled_old:
            if LEDamount > 0:
                pmgain[jled] = dxCh*np.float32(LEDsumm)/np.float32(LEDamount)
                if LEDsumm > 0:
                    LEDcoeff = np.float32(led_ref)/pmgain[jled]
            pulseheight[jtmark: jpul] *= LEDcoeff
            jtmark = jpul
            LEDsumm = 0
            LEDamount = 0
            LEDcoeff = 0
        if flg_led[jpul]: # LED single pulse
            LEDsumm += totalintegral[jpul]
            LEDamount += 1
        jled_old = jled

    return pmgain, pulseheight

@nb.njit
def BaselineCond2(bl_start, max_diff, pulses, pulse_len, maxpos, max_lg):

    n_pulses = maxpos.shape[0]
    pulse_basestart = pulse_len - bl_start
    blstart_h = bl_start//2
    totalintegral = np.zeros(n_pulses, dtype=np.float32)

    for jpul in range(n_pulses):
        newpulse_len = 0
        pulse = pulses[jpul]
        if pulse_basestart[jpul] >= maxpos[jpul]:
            newpulse_len = pulse_len[jpul] - blstart_h
        else:
            aver1 = 0.
            for i in range(bl_start):
                aver1 += pulse[i]
            aver1 /= float(bl_start)
            for j in range(maxpos[jpul], pulse_basestart[jpul]):
                aver2 = 0.
                for i in range(j, j + bl_start):
                    aver2 += pulse[i]
                aver2 /= float(bl_start)
                if np.abs(aver2 - aver1) < max_diff:
                    if max_lg[jpul] > j + blstart_h:
                        newpulse_len = max_lg[jpul]
                    else:
                        newpulse_len = j + blstart_h
                    break
                if (j == pulse_basestart[jpul] - 1) :
                    newpulse_len = pulse_len[jpul] - blstart_h

        for j in range(1, newpulse_len-1):
            totalintegral[jpul] += pulse[j]
        totalintegral[jpul] += 0.5*(pulse[0] + pulse[newpulse_len-1])

    return totalintegral


@nb.njit
def Baseline(basestart, baseend, pulse_len, pulses):

    n_pulses = pulses.shape[0]
    pulse_baseend = pulse_len - baseend
    baseline = np.zeros(n_pulses, dtype=np.float32)
    for jpul in range(n_pulses):
# Baseline subtraction
        nind = 0
        for j in range(basestart):
            baseline[jpul] += pulses[jpul, j]
            nind += 1
        for j in range(pulse_baseend[jpul], pulse_len[jpul]):
            if j >= basestart:
                baseline[jpul] += pulses[jpul, j]
                nind += 1
        baseline[jpul] /= float(nind)
    return baseline

@nb.njit
def PileUpDet(nfront, ntail, nthres, front_led, tail_led, flags, pulses):

    n_pulses, pulse_len = pulses.shape
    flg_peaks = np.zeros(n_pulses, dtype=np.int32)

    for jpul, flg in enumerate(flags):
        pulse = pulses[jpul]
        if flg:
            pulse_width = front_led + tail_led
            pulse_front = pulse[front_led : -tail_led] - nthres
        else:
            pulse_width = nfront + ntail
            pulse_front = pulse[nfront : -ntail] - nthres

        pulse_max = np.maximum(pulse[:-pulse_width], pulse[pulse_width:])
        flg_peaks[jpul] = 0
        jt = 0
        while jt < pulse_len-pulse_width:
            if pulse_front[jt] > pulse_max[jt]:
                flg_peaks[jpul] += 1
                jt += pulse_width
            jt += 1
        flg_peaks[jpul] 

    return flg_peaks


class DPSD:


    def __init__(self, dic_in, t_ranges=None):

        self.d_int = {}
        self.d_flt = {}
        int_list = ('BaselineStart', 'BaselineStart', 'BaselineEnd', \
            'DDlower', 'DDupper', 'DTlower', 'DTupper', \
            'ShortGate', 'LongGate', 'Threshold', 'Front', 'Tail', \
            'Marker', 'SaturationHigh', 'SaturationLow', 'ToFWindowLength', \
            'PH_nChannels', 'PS_nChannels', 'LineChange', \
            'LEDxmin', 'LEDxmax', 'LEDymin', 'LEDymax', 'LEDFront', 'LEDTail', 'LEDreference')
        flt_list = ('LEDdt', 'TimeBin', 'TBeg', 'TEnd', 'MaxDifference', 'Slope1', 'Slope2', 'Offset', )

        for key in int_list:
            self.d_int[key] = int(dic_in[key])
        for key in flt_list:
            self.d_flt[key] = float(dic_in[key])
        HAfile = dic_in['HAfile'].strip()
        if HAfile != '':
            self.HAfile = HAfile
            self.run(HAfile, t_ranges=t_ranges)
        else:
            n_shots = np.atleast_1d(eval(str(dic_in['Shots'])))
            for nshot in n_shots:
                logger.info(nshot)
                self.nshot = int(nshot)
                shot100 = self.nshot//100
                filepath = '/afs/ipp/augd/rawfiles/NSP/%d/%d' %(shot100, self.nshot)
                HAfile = '%s/HA_%d.dat' %(filepath, self.nshot)
                self.HAfile = HAfile
                self.run(HAfile, t_ranges=t_ranges, check_md5=True)
                if 'SFwrite' in dic_in.keys():
                    if dic_in['SFwrite']:
                        self.sfwrite(exp=dic_in['SFexp'])


    def run(self, HAfile, t_ranges=None, check_md5=False):

        nxCh = self.d_int['PH_nChannels']
        nyCh = self.d_int['PS_nChannels']
        dxCh = np.float32(nxCh)/np.float32(self.d_int['Marker'])

        min_winlen = max(self.d_int['BaselineStart'], self.d_int['BaselineEnd'])
        ha = read_ha.READ_HA(HAfile, check_md5=check_md5, min_winlen=min_winlen, max_winlen=self.d_int['ToFWindowLength'])
        if not ha.status:
            return

        if t_ranges is None:
            if self.d_flt['TEnd'] <= 0:
                self.d_flt['TEnd'] = ha.t_events[-1] # Take all time events
            (tind, ) = np.where((ha.t_events >= self.d_flt['TBeg']) & (ha.t_events <= self.d_flt['TEnd']))
            self.dt = self.d_flt['TEnd'] - self.d_flt['TBeg']
        else: # Force time ranges (make sure they don't overlap!)
            depth = lambda L: isinstance(L, list) and max(map(depth, L))+1 # list depth
            if depth(t_ranges) == 1:
                t_ranges = [t_ranges]    
            tind = np.array([], dtype=np.int32)
            self.dt = 0
            for jint, t_ran in enumerate(t_ranges):
                (ind, ) = np.where((ha.t_events >= t_ran[0]) & (ha.t_events <= t_ran[1]))
                tind = np.append(tind, ind)
                self.dt += t_ran[1] - t_ran[0]

        self.time = ha.t_events[tind]
        logger.info('TBeg = %8.4f' %self.time[0]) 
        logger.info('TEnd = %8.4f' %self.time[-1]) 
        n_pulses = len(self.time)
        n_timebins = int((self.time[-1] - self.time[0])/self.d_flt['TimeBin'])
        n_led = int((self.time[-1] - self.time[0])/self.d_flt['LEDdt'])

        self.time_cnt = self.time[0] + self.d_flt['TimeBin']*(0.5 + np.arange(n_timebins))
        self.time_led = self.time[0] + self.d_flt['LEDdt']*(0.5 + np.arange(n_led))

        self.winlen = ha.winlen[tind]
        pulses = ha.pulses[tind]

# Initialise
        flg_sat   = np.zeros(n_pulses)
        self.flg_peaks = np.zeros(n_pulses)

# =========== 
# Time loop
# =========== 

        logger.info('Starting time loop')

        logger.info('# pulses: %d' %n_pulses)

        maxpos = np.argmax(pulses, axis=1)
        pulse_max = np.max(pulses, axis=1)

        max_LG = np.minimum(maxpos + self.d_int['LongGate' ], self.winlen)
        max_SG = np.minimum(maxpos + self.d_int['ShortGate'], self.winlen)
        tof_win_len = self.d_int['ToFWindowLength']
        sat_high = float(self.d_int['SaturationHigh'])
        sat_low  = float(self.d_int['SaturationLow'])
        pulse_len = np.minimum(self.winlen, tof_win_len)
        pulse_baseend = pulse_len - self.d_int['BaselineEnd']
        
        logger.info('Baseline subtraction')
        self.pulses = pulses.astype(np.float32)
        baseline = Baseline(self.d_int['BaselineStart'], self.d_int['BaselineEnd'], pulse_len, self.pulses)
        self.pulses -= baseline[:, None]

# Saturation detection
        logger.info('Saturation detection')
        (ind_sat_high, ) = np.where(np.max(self.pulses, axis=1) > sat_high)
        (ind_sat_low, )  = np.where(np.min(self.pulses, axis=1) < sat_low) 
        flg_sat[ind_sat_high] = 1
        flg_sat[ind_sat_low]  = 2

        logger.info('Baseline conditioned 2') 

        self.TotalIntegral = BaselineCond2(self.d_int['BaselineStart'], self.d_flt['MaxDifference'], self.pulses, pulse_len, maxpos, max_LG)
        self.ShortIntegral = slice_trapz(self.pulses, maxpos, max_SG)
        self.LongIntegral  = slice_trapz(self.pulses, maxpos, max_LG)
        ind3 = np.where(self.LongIntegral > 0)[0]

        self.PulseHeight = dxCh*self.TotalIntegral
        self.PulseShape = np.zeros(n_pulses, dtype=np.float32)
        self.PulseShape[ind3] = np.float32(nyCh)*self.ShortIntegral[ind3]/self.LongIntegral[ind3]

# LED evaluation

        self.flg = {}

        self.flg['led'] =  \
            (self.PulseHeight > float(self.d_int['LEDxmin'])) & \
            (self.PulseHeight < float(self.d_int['LEDxmax'])) & \
            (self.PulseShape > float(self.d_int['LEDymin'])) & \
            (self.PulseShape < float(self.d_int['LEDymax'])) 

        logger.info('Pile-up detection')

        self.flg_peaks = PileUpDet(self.d_int['Front'], self.d_int['Tail'], self.d_int['Threshold'], self.d_int['LEDFront'], self.d_int['LEDTail'], self.flg['led'], pulses)

# LED correction

        logger.info('LED correction')

        self.pmgain, self.PulseHeight = led_correction(self.d_flt['LEDdt'], dxCh, self.d_int['LEDreference'], self.time, self.TotalIntegral, self.flg['led'])

        self.TotalIntegral = self.PulseHeight/dxCh

        flg_slope1 = (self.PulseHeight <= self.d_int['LineChange'])
        flg1n = (self.PulseShape <= self.d_flt['Offset'] + self.d_flt['Slope1']*self.PulseHeight)
        offset2 = self.d_flt['Offset'] + self.d_flt['Slope1']*self.d_int['LineChange']
        flg2n = (self.PulseShape <= offset2 + self.d_flt['Slope2']*(self.PulseHeight-self.d_int['LineChange']))

        flg1  =   flg_slope1  & flg1n
        flg2  = (~flg_slope1) & flg2n
        flg1g =   flg_slope1  & (~flg1n)
        flg2g = (~flg_slope1) & (~flg2n)
        self.flg['sat'] = (flg_sat > 0)
        self.flg['pileup'] = (self.flg_peaks > 1)
        self.flg['phys'] =  (~self.flg['sat']) & (~self.flg['led']) & (~self.flg['pileup'])
        self.flg['neut1']  = (flg1  + flg2 ) & (self.flg['phys'])
        self.flg['gamma1'] = (flg1g + flg2g) & (self.flg['phys'])
        self.flg['DD'] = self.flg['neut1'] & \
            (self.PulseHeight >= self.d_int['DDlower']) & \
            (self.PulseHeight <= self.d_int['DDupper'])
        self.flg['DT'] = self.flg['neut1'] & \
            (self.PulseHeight >= self.d_int['DTlower']) & \
            (self.PulseHeight <= self.d_int['DTupper'])

        self.event_type = np.zeros(n_pulses, dtype=np.int32) - 1
        self.event_type[self.flg['neut1']]  = 0
        self.event_type[self.flg['gamma1']] = 1
        self.event_type[self.flg['pileup']] = 2
        self.event_type[self.flg['led']]    = 3

        cnt_list = ('neut1', 'gamma1', 'led', 'pileup', 'sat', 'phys', 'DD', 'DT')
        nxCh = self.d_int['PH_nChannels']

        self.cnt = {}
        self.phs = {}

        for spec in cnt_list:
            cnt, tedges = np.histogram(self.time[self.flg[spec]], bins=n_timebins, range=[self.time_cnt[0]-0.5*self.d_flt['TimeBin'], self.time_cnt[-1]+0.5*self.d_flt['TimeBin']])
            phs, edges = np.histogram(np.float32(nxCh)/np.float32(self.d_int['Marker'])*self.TotalIntegral[self.flg[spec]], bins=nxCh, range=[-0.5, nxCh + 0.5])
            self.cnt[spec] = cnt.astype(np.float32)
            self.phs[spec] = phs.astype(np.float32)/self.dt
            logger.info('%s %d', spec, np.sum(self.flg[spec]))

# Move to 1/s units
        for spec in self.cnt.keys():
            self.cnt[spec] /= self.d_flt['TimeBin']

        total = self.cnt['neut1'] + self.cnt['gamma1'] + self.cnt['led']
# Assuming pile-ups are all 2 events per window
        pup_frac = 1 + 2.*self.cnt['pileup']/total
        self.cnt['neut2' ] = pup_frac*self.cnt['neut1'] 
        self.cnt['gamma2'] = pup_frac*self.cnt['gamma1']


    def sfwrite(self, fsfh='%s/NSP00000.sfh' %dpsd_dir, exp='AUGD'):

        import aug_sfutils as sf
        from aug_sfutils import sfhmod
        ww = sf.WW()

        diag = 'NSP'
        nsp = sf.SFREAD('NSP', self.nshot, exp=exp)
        if nsp.status:
            logger.error('NSP shotfile for #%d exists already' %self.nshot)
            return

        nt = len(self.time_cnt)
        sfh = sfhmod.SFHMOD(fin=fsfh)
        for lbl in ['time'] + sig1d:
            sfh.modtime(lbl, nt)
        sfh.write(fout=fsfh)

        os.chdir(dpsd_dir)
        if ww.Open(exp, diag, self.nshot):
            status = ww.SetSignal('time', np.array(self.time_cnt, dtype=np.float32))
            for lbl in sig1d:
                status = ww.SetSignal(lbl, np.array(self.cnt[lbl], dtype=np.float32))
            ww.Close()
