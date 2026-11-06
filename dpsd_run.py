import sys, os, logging, time, shutil

import numpy as np
import numba as nb
import read_ha
from pathlib import Path

logger = logging.getLogger('dpsd.run')

sig1d = ['time', 'neut1', 'neut2', 'gamma1', 'gamma2', 'led', 'pileup']
dpsd_dir = Path(__file__).resolve().parent

@nb.njit
def slice_trapz(a, bnd_l, bnd_r):
    b = np.empty(a.shape[0])
    for j in range(a.shape[0]):
        b[j] = np.sum(a[j, bnd_l[j]+1: bnd_r[j]-1])
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
            aver1 = np.mean(pulse[:bl_start])
            for j in range(maxpos[jpul], pulse_basestart[jpul]):
                aver2 = np.mean(pulse[j: j+bl_start])
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
    baseline = np.sum(pulses[:, :basestart], axis=1)
    for jpul in range(n_pulses):
        nind = basestart
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
    return flg_peaks


class DPSD:


    def __init__(self, dic_in, t_ranges=None, nesp_root=Path('/shares/experiments/aug-rawfiles/NSP')):

        self.status = True

        self.setup = dic_in
        io_d = dic_in['io']
        HAfile = io_d['HA*.dat file'].strip()
        if HAfile != '':
            self.run(Path(HAfile), t_ranges=t_ranges)
        else:
            n_shots = np.atleast_1d(eval(str(io_d['Shots'])))
            for nshot in n_shots:
                logger.info(nshot)
                self.nshot = int(nshot)
                sshot = f'{nshot:05d}'
                shot100 = sshot[:3]
                HAfile = nesp_root / shot100 / sshot / f'HA_{sshot}.dat'
                if not HAfile.is_file():
                    HAfile = nesp_root / shot100 / f'HA_{sshot}.dat'
                self.run(HAfile, t_ranges=t_ranges, check_md5=True)
                if 'Write shotfiles' in io_d.keys():
                    if self.status and io_d['Write shotfiles']:
                        self.sfwrite(exp=io_d['Shotfile exp'], force=io_d['Force SF write'])


    def run(self, HAfile, t_ranges=None, check_md5=False):

        nxCh = self.setup['separation']['#bins Pulse Height']
        nyCh = self.setup['separation']['#bins Pulse Shape']
        dxCh = np.float32(nxCh)/np.float32(self.setup['separation']['Marker'])

        min_winlen = max(self.setup['peak']['Baseline start'], self.setup['peak']['Baseline end'])
        print(type(min_winlen), type(self.setup['setup']['#samples for analysis']))
        ha = read_ha.READ_HA(HAfile, check_md5=check_md5, min_winlen=min_winlen, max_winlen=self.setup['setup']['#samples for analysis'])
        self.status = ha.status
        if not self.status:
            return

        if t_ranges is None:
            if self.setup['setup']['End time'] <= 0:
                self.setup['setup']['End time'] = ha.t_events[-1] # Take all time events
            (tind, ) = np.where((ha.t_events >= self.setup['setup']['Start time']) & (ha.t_events <= self.setup['setup']['End time']))
            self.dt = self.setup['setup']['End time'] - self.setup['setup']['Start time']
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
        logger.info('Start time = %8.4f', self.time[0]) 
        logger.info('End time = %8.4f', self.time[-1]) 
        n_pulses = len(self.time)
        n_timebins = int((self.time[-1] - self.time[0])/self.setup['setup']['Time step'])
        n_led = int((self.time[-1] - self.time[0])/self.setup['led']['LED time sampling'])

        self.cnt = {'time': self.time[0] + self.setup['setup']['Time step']*(0.5 + np.arange(n_timebins))}
        self.time_led = self.time[0] + self.setup['led']['LED time sampling']*(0.5 + np.arange(n_led))

        self.winlen = ha.winlen[tind]
        pulses = ha.pulses[tind]

# Initialise
        self.flg_peaks = np.zeros(n_pulses)

# =========== 
# Time loop
# =========== 

        logger.info('Starting time loop')

        logger.info('# pulses: %d', n_pulses)

        maxpos = np.argmax(pulses, axis=1)
        pulse_max = np.max(pulses, axis=1)

        max_LG = np.minimum(maxpos + self.setup['peak']['Long gate' ], self.winlen)
        max_SG = np.minimum(maxpos + self.setup['peak']['Short gate'], self.winlen)
        tof_win_len = self.setup['setup']['#samples for analysis']
        sat_high = float(self.setup['peak']['Saturation upper limit'])
        sat_low  = float(self.setup['peak']['Saturation lower limit'])
        pulse_len = np.minimum(self.winlen, tof_win_len)
        pulse_baseend = pulse_len - self.setup['peak']['Baseline end']

        logger.info('Baseline subtraction')
        self.pulses = pulses.astype(np.float32)
        baseline = Baseline(self.setup['peak']['Baseline start'], self.setup['peak']['Baseline end'], pulse_len, self.pulses)
        self.pulses -= baseline[:, None]

# Saturation detection
        logger.info('Saturation detection')
        (ind_sat_high, ) = np.where(np.max(self.pulses, axis=1) > sat_high)
        (ind_sat_low, )  = np.where(np.min(self.pulses, axis=1) < sat_low )
        flg_sat   = np.zeros(n_pulses)
        flg_sat[ind_sat_high] = 1
        flg_sat[ind_sat_low]  = 2

        logger.info('Baseline conditioned 2') 

        self.TotalIntegral = BaselineCond2(self.setup['peak']['Baseline start'], self.setup['peak']['Maximum difference'], self.pulses, pulse_len, maxpos, max_LG)
        self.ShortIntegral = slice_trapz(self.pulses, maxpos, max_SG)
        self.LongIntegral  = slice_trapz(self.pulses, maxpos, max_LG)
        ind3 = np.where(self.LongIntegral > 0)[0]

        self.PulseHeight = dxCh*self.TotalIntegral
        self.PulseShape = np.zeros(n_pulses, dtype=np.float32)
        self.PulseShape[ind3] = np.float32(nyCh)*self.ShortIntegral[ind3]/self.LongIntegral[ind3]

# LED evaluation

        self.flg = {}

        self.flg['led'] =  \
            (self.PulseHeight > float(self.setup['led']['Min PH bin for LED detection'])) & \
            (self.PulseHeight < float(self.setup['led']['Max PH bin for LED detection'])) & \
            (self.PulseShape  > float(self.setup['led']['Min PS bin for LED detection'])) & \
            (self.PulseShape  < float(self.setup['led']['Max PS bin for LED detection']))

        logger.info('Pile-up detection')

        self.flg_peaks = PileUpDet(self.setup['peak']['Front'], self.setup['peak']['Tail'], self.setup['peak']['Threshold'], self.setup['led']['LED front'], self.setup['led']['LED tail'], self.flg['led'], pulses)

# LED correction

        logger.info('LED correction')

        self.pmgain, self.PulseHeight = led_correction(self.setup['led']['LED time sampling'], dxCh, self.setup['led']['LED reference bin'], self.time, self.TotalIntegral, self.flg['led'])

        self.TotalIntegral = self.PulseHeight/dxCh

        flg_slope1 = (self.PulseHeight <= self.setup['separation']['Bin line1 -> line2'])
        flg1n = (self.PulseShape <= self.setup['separation']['Offset of 1st sep.line'] + self.setup['separation']['Slope of 1st sep.line']*self.PulseHeight)
        offset2 = self.setup['separation']['Offset of 1st sep.line'] + self.setup['separation']['Slope of 1st sep.line']*self.setup['separation']['Bin line1 -> line2']
        flg2n = (self.PulseShape <= offset2 + self.setup['separation']['Slope of 2nd sep.line']*(self.PulseHeight-self.setup['separation']['Bin line1 -> line2']))

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
            (self.PulseHeight >= self.setup['separation']['Lower PH-limit for DD']) & \
            (self.PulseHeight <= self.setup['separation']['Upper PH-limit for DD'])
        self.flg['DT'] = self.flg['neut1'] & \
            (self.PulseHeight >= self.setup['separation']['Lower PH-limit for DT']) & \
            (self.PulseHeight <= self.setup['separation']['Upper PH-limit for DT'])

        self.event_type = np.zeros(n_pulses, dtype=np.int32) - 1
        self.event_type[self.flg['neut1']]  = 0
        self.event_type[self.flg['gamma1']] = 1
        self.event_type[self.flg['pileup']] = 2
        self.event_type[self.flg['led']]    = 3

        cnt_list = ('neut1', 'gamma1', 'led', 'pileup', 'sat', 'phys', 'DD', 'DT')
        nxCh = self.setup['separation']['#bins Pulse Height']

        self.phs = {}

        for spec in cnt_list:
            cnt, tedges = np.histogram(self.time[self.flg[spec]], bins=n_timebins, range=[self.cnt['time'][0]-0.5*self.setup['setup']['Time step'], self.cnt['time'][-1]+0.5*self.setup['setup']['Time step']])
            phs, edges = np.histogram(np.float32(nxCh)/np.float32(self.setup['separation']['Marker'])*self.TotalIntegral[self.flg[spec]], bins=nxCh, range=[-0.5, nxCh + 0.5])
            self.cnt[spec] = cnt.astype(np.float32)
            self.phs[spec] = phs.astype(np.float32)/self.dt
            logger.info('%s %d', spec, np.sum(self.flg[spec]))

# Move to 1/s units
        for spec in self.cnt.keys():
            self.cnt[spec] /= self.setup['setup']['Time step']

        total = self.cnt['neut1'] + self.cnt['gamma1'] + self.cnt['led']
# Assuming pile-ups are all 2 events per window
        pup_frac = 1 + 2.*self.cnt['pileup']/total
        self.cnt['neut2' ] = pup_frac*self.cnt['neut1'] 
        self.cnt['gamma2'] = pup_frac*self.cnt['gamma1']


    def sfwrite(self, exp='augd', force=False):

        import aug_sfutils as sf
        from aug_sfutils import sfhmod

        diag = 'NSP'
        nsp = sf.SFREAD('NSP', self.nshot, exp=exp)
        if nsp.status and not force:
            logger.info('NSP shotfile for #%d exists already' %self.nshot)
            return

        counts = {key: self.cnt[key] for key in sig1d}
        counts['neutD-T'] = self.cnt['DT']
        shutil.copy(dpsd_dir / 'NSP00000.sfh.temp', dpsd_dir / 'NSP00000.sfh')
        status = sf.write_sf(self.nshot, counts, dpsd_dir, 'NSP', exp=exp)
