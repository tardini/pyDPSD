import logging
import numpy as np


fmt = logging.Formatter('%(asctime)s | %(name)s | %(levelname)s: %(message)s', '%H:%M:%S')
logger = logging.getLogger('read_HA')
hnd = logging.StreamHandler()
hnd.setFormatter(fmt)
logger.addHandler(hnd)
logger.setLevel(logging.DEBUG)


class READ_HA:


    def __init__(self, fin, min_winlen=0):


        logger.info('Reading binary %s', fin)
        data = np.fromfile(fin, dtype=np.uint16)
        tdiff = []
        self.boundaries = []
        skipped_odd = []
        (ind, ) = np.where((data >= 0) & (data < 3))
        logger.info('Getting t_diff and win_len')

        for jpos in ind:
            if (data[jpos+2] == data[jpos]) and (data[jpos+1] == data[jpos+3]-1):
                self.boundaries.append(jpos)
                if data[jpos] == 0:
                    tdiff.append(data[jpos+3])
                elif data[jpos] == 1:
                    tdiff.append(32768 + data[jpos+3])
                elif data[jpos] == 2:
                    tdiff.append(65536 + data[jpos+3])
        self.boundaries.append(len(data)) # Retain final pulse too, unlike *.bin
        winlen = np.diff(self.boundaries) - 4

        (ind_odd, ) = np.where(winlen %2 == 1)
        (ind_wneg, ) = np.where(winlen < 0)
        logger.debug('Skipped %d pulses with odd window length' , len(ind_odd))
        logger.info('Skipped %d pulses with window length <= 0', len(ind_wneg))

        (ind_ok, ) = np.where((winlen %2 == 0) & (winlen > min_winlen))
        self.winlen = winlen[ind_ok]
        self.tdiff = np.array(tdiff)[ind_ok]
        win_start = np.array(self.boundaries[:-1])[ind_ok]

        rawdata = data.astype(np.int16)
        rawdata -= 32768
        (ind_neg, ) = np.where(rawdata > 8192)
        rawdata[ind_neg] -= 16384
        rawdata *= -1
# Entry-inversion observed by Luca Giacomelli
        ndat = len(rawdata)
        self.rawdata = np.append(rawdata[1::2], rawdata[::2]).reshape((2,  ndat//2)).T.ravel()

        self.pulses   = []
        logger.info('Reading pulses')
        for jwin, jpos in enumerate(win_start):
            self.pulses.append(self.rawdata[jpos + 4: jpos + 4 + self.winlen[jwin]])

        logger.debug('Min winlen %d %d', np.min(winlen), np.min(self.winlen)) 
        logger.debug('%d', len(self.pulses))
        self.t_events = 1e-8*(np.cumsum(self.tdiff, dtype=np.float32))


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
