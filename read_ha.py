import os, logging
import numpy as np
import numba as nb

fmt = logging.Formatter('%(asctime)s | %(name)s | %(levelname)s: %(message)s', '%H:%M:%S')
logger = logging.getLogger('read_HA')
hnd = logging.StreamHandler()
hnd.setFormatter(fmt)
logger.addHandler(hnd)
logger.setLevel(logging.INFO)


@nb.njit
def raw2pulse(max_winlen, win_start, pulse_len, rawdata):

    win_end = win_start + pulse_len
    pulses = np.zeros((win_start.shape[0], max_winlen))

# Determine shift using a minimum-derivative**2 approach, pulse by pulse

    for jwin, jpos in enumerate(win_start):
        pulse = rawdata[jpos : win_end[jwin]]
        pulse_even = pulse[:-1:2]
        pulse_odd  = pulse[1::2]
        len0 = len(pulse_odd)
        min_tens = 1e8
        for j in range(1, 3):
            pulse2 = np.stack((pulse_even[j:], pulse_odd[:len0-j])).T.ravel()
            der2 = (pulse2[1:] - pulse2[:-1])**2
            tension = 0.
            for x in der2:
                tension += x
            if tension < min_tens:
                pulse_ok = pulse2
                min_tens = tension

        for j in range(3):
            pulse2 = np.stack((pulse_odd[j:], pulse_even[:len0-j])).T.ravel()
            der2 = (pulse2[1:] - pulse2[:-1])**2
            tension = 0.
            for x in der2:
                tension += x
            if tension < min_tens:
                pulse_ok = pulse2
                min_tens = tension

        len_pul = len(pulse_ok)
        pulses[jwin][: len_pul] = pulse_ok

    return pulses


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
        (boundaries, ) = np.where( (data[:-3] >= 0) & (data[:-3] < 3) & \
            (data[2:-1] < 3) & (data[2:-1] >= 0) & (data1[1:-2] == data[3:]) )

        data32 = data.astype(np.uint32)
        tdiff = data32[boundaries + 3] + data32[boundaries]*32768

        self.boundaries = np.append(boundaries, len(data)) # Retain final pulse too, unlike *.bin
        winlen = np.diff(self.boundaries) - 4

        (ind_odd, ) = np.where(winlen %2 == 1)
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
# Entry-inversion observed by Luca Giacomelli
        n_pulses = len(self.winlen)


        if max_winlen is None:
            max_winlen = np.max(self.winlen)

        logger.info('Reading pulses')
        pulse_len = np.minimum(self.winlen, max_winlen)

        self.pulses = raw2pulse(max_winlen, win_start, pulse_len, data)

        self.t_events = 1e-8*(np.cumsum(tdiff, dtype=np.float32))[ind_ok]
        logger.debug('Min winlen %d %d', np.min(winlen), np.min(self.winlen)) 
        logger.debug('%d', len(self.pulses))
