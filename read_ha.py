import os, logging
import numpy as np
import numba as nb
from numba.typed import List

fmt = logging.Formatter('%(asctime)s | %(name)s | %(levelname)s: %(message)s', '%H:%M:%S')
logger = logging.getLogger('read_HA')
hnd = logging.StreamHandler()
hnd.setFormatter(fmt)
logger.addHandler(hnd)
logger.setLevel(logging.INFO)


@nb.njit
def minTension(pulse_in):
    pulse_even = np.ascontiguousarray(pulse_in[:-1:2])
    pulse_odd  = np.ascontiguousarray(pulse_in[1::2])
    len0 = len(pulse_odd)
    min_tens = 1e8
    pulse_ok = np.empty(0, dtype=pulse_in.dtype)
    jmin = -1
    for j in range(3):
        N = len0 - j
        pulse2 = np.empty(2 * N, dtype=pulse_in.dtype)
        pulse2[0::2] = pulse_odd[j:]
        pulse2[1::2] = pulse_even[:N]
        tension = np.sum((pulse2[1:] - pulse2[:-1])**2)
        if tension < min_tens:
            pulse_ok = pulse2
            min_tens = tension
            jmin = j
    return jmin, pulse_ok

@nb.njit
def raw2pulse(max_winlen, win_start, pulse_len, rawdata):
    win_end = win_start + pulse_len
    n_pulses = win_start.shape[0]
    pulses = np.zeros((n_pulses, max_winlen))
    indBad = List()
# Determine shift using a minimum-derivative**2 approach, pulse by pulse
    for jwin in range(n_pulses):
        jpos = win_start[jwin]
        pulse = rawdata[jpos : win_end[jwin]]
        jmin, pulse_ok = minTension(pulse)
        if jmin > 0:
            indBad.append(jwin)
        len_pul = len(pulse_ok)
        pulses[jwin, : len_pul] = pulse_ok
    return indBad, pulses


class READ_HA:


    def __init__(self, fin, check_md5=False, min_winlen=0, max_winlen=None):

        self.status = True
        logger.info('Reading binary %s', fin)
        fmd5 = fin + '.md5'
        if not os.path.isfile(fin):
            logger.error('File %s not found', fin)
            self.status = False
            return
        if check_md5:
            if not os.path.isfile(fmd5): # ensures integrity of fin
                logger.error('File %s not found', fmd5)
                self.status = False
                return
        data = np.fromfile(fin, dtype=np.uint16)

        logger.info('Getting t_diff and win_len')
        data1 = data + 1
        (boundaries, ) = np.where(
            np.isin(data[ :-3], [0, 1, 2]) & \
            np.isin(data[2:-1], [0, 1, 2]) & \
            (data1[1:-2] == data[3:]) )

        data32 = data.astype(np.uint32)
        tdiff = data32[boundaries + 3] + data32[boundaries]*32768

        self.boundaries = np.append(boundaries, len(data)) # Retain final pulse too, unlike *.bin
        winlen = np.diff(self.boundaries) - 4

        (ind_odd , ) = np.where(winlen %2 == 1)
        (ind_wneg, ) = np.where(winlen < 0)
        logger.debug('Skipped %d pulses with odd window length' , len(ind_odd))
        logger.info('Skipped %d pulses with window length <= 0', len(ind_wneg))

        (ind_ok, ) = np.where((winlen %2 == 0) & (winlen > min_winlen))
        self.winlen = winlen[ind_ok]

        win_start = boundaries[ind_ok] + 4

        data = data.astype(np.int16)
        data -= 32768
        (ind_neg, ) = np.where(data > 8192)
        data[ind_neg] -= 16384
        data *= -1

        n_pulses = len(self.winlen)

        if max_winlen is None:
            max_winlen = np.max(self.winlen)

        pulse_len = np.minimum(self.winlen, max_winlen)

# Entry-inversion observed by Luca Giacomelli
        logger.info('Sorting faulty ADC synchronisation')
        indBad, self.pulses = raw2pulse(max_winlen, win_start, pulse_len, data)
        n_sorted = len(indBad)
        logger.info('Sorted ADC for %d points out of %d', n_sorted, win_start.shape[0])

        self.t_events = 1e-8*(np.cumsum(tdiff, dtype=np.float32))[ind_ok]
        logger.debug('Min winlen %d %d', np.min(winlen), np.min(self.winlen)) 
        logger.debug('%d', len(self.pulses))

#        import matplotlib.pylab as plt
#        jbad = indBad[0]
#        plt.plot(data[win_start[jbad]: win_start[jbad] + pulse_len[jbad]])
#        plt.plot(self.pulses[jbad, :])
#        plt.show()
