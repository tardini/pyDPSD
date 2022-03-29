import logging
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
    for jwin, jpos in enumerate(win_start):
        pulses[jwin][: pulse_len[jwin]] = rawdata[jpos : win_end[jwin]]
    return pulses


class READ_HA:


    def __init__(self, fin, min_winlen=0, max_winlen=None):


        logger.info('Reading binary %s', fin)
        data = np.fromfile(fin, dtype=np.uint16)

        logger.info('Getting t_diff and win_len')
        data1 = data + 1
        (boundaries, ) = np.where( (data[:-3] >= 0) & (data[:-3] < 3) & \
            (data[2:-1] < 3) & (data[2:-1] >= 0) & (data1[1:-2] == data[3:]) )

        tdiff = data[boundaries + 3] + data[boundaries]*32768

        self.boundaries = np.append(boundaries, len(data)) # Retain final pulse too, unlike *.bin
        winlen = np.diff(self.boundaries) - 4

        (ind_odd, ) = np.where(winlen %2 == 1)
        (ind_wneg, ) = np.where(winlen < 0)
        logger.debug('Skipped %d pulses with odd window length' , len(ind_odd))
        logger.info('Skipped %d pulses with window length <= 0', len(ind_wneg))

        (ind_ok, ) = np.where((winlen %2 == 0) & (winlen > min_winlen))
        self.winlen = winlen[ind_ok]
        self.tdiff = tdiff[ind_ok]
        win_start = boundaries[ind_ok] + 4

        data = data.astype(np.int16)
        data -= 32768
        (ind_neg, ) = np.where(data > 8192)
        data[ind_neg] -= 16384
        data *= -1
# Entry-inversion observed by Luca Giacomelli
        n_pulses = len(self.winlen)
        self.rawdata = np.stack((data[1::2], data[::2])).T.ravel()

        if max_winlen is None:
            max_winlen = np.max(self.winlen)

        logger.info('Reading pulses')
        pulse_len = np.minimum(self.winlen, max_winlen)

        self.pulses = raw2pulse(max_winlen, win_start, pulse_len, self.rawdata)

        self.t_events = 1e-8*(np.cumsum(self.tdiff, dtype=np.float32))
        logger.debug('Min winlen %d %d', np.min(winlen), np.min(self.winlen)) 
        logger.debug('%d', len(self.pulses))


if __name__ == '__main__':

    import matplotlib.pylab as plt

    nshot = 40490
    nshot = 40525

    fin = 'HA_%d.dat' %nshot
    ha = READ_HA(fin)

    tfile = '%dt.bin' %nshot
    lfile = '%dl.bin' %nshot
    dfile = '%d.bin'  %nshot

    logger.info('Reading bin files')
    tdiff  = np.fromfile(tfile, dtype = np.uint32)
    winlen = np.fromfile(lfile, dtype = np.uint16)
    data   = np.fromfile(dfile, dtype = np.int16)
    logger.info('Done bin files')

    print('Tdiff')
    print(len(ha.tdiff), len(tdiff))
    if len(ha.tdiff) == len(tdiff):
        print(np.max(np.abs(ha.tdiff - tdiff)))

    print('Winlen')
    print(len(ha.winlen), len(winlen))
    if len(ha.winlen) == len(winlen) + 1:
        print(np.max(np.abs(ha.winlen[:-1] - winlen)))

    print('Data')

    pos = 0

    pulse_ha_test1 = []
    pulse_ha_test2 = []
    pulse_test1 = []
    pulse_test2 = []
    for jpul, pulseha in enumerate(ha.pulses):
        if jpul == len(winlen):
            break
        pulse = data[pos: pos + winlen[jpul]]
        pos += winlen[jpul]
        if len(pulseha) != len(pulse):
            pulse_ha_test1.append(pulseha)
            pulse_test1.append(pulse)
        elif np.max(np.abs(pulseha - pulse)) > 0:
            pulse_ha_test2.append(pulseha)
            pulse_test2.append(pulse)

    print(len(pulse_ha_test1), len(pulse_ha_test2), len(ha.pulses))
    nplot = 50

    plt.figure(1, (14., 6.4))
    plt.subplot(2, 2, 1)
    for pulseha in pulse_ha_test1[:nplot]:
        plt.plot(pulseha)
    plt.subplot(2, 2, 2)
    for pulse in pulse_test1[:nplot]:
        plt.plot(pulse)
    plt.subplot(2, 2, 3)
    for pulseha in pulse_ha_test2[:nplot]:
        plt.plot(pulseha)
    plt.subplot(2, 2, 4)
    for pulse in pulse_test2[:nplot]:
        plt.plot(pulse)
    plt.show()

    
#    data = np.fromfile(fin, dtype=np.uint16)
#    print(data[    : 1000])
#    print(data[1000: 2000])
#    print(data[2000: 3000])
