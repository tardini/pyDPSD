import os, json
import dpsd_run
import numpy as np
import matplotlib.pylab as plt

locdir = os.path.dirname(os.path.realpath(__file__))
f_json = '%s/settings/moko.json' %locdir
with open(f_json, 'r') as fjson:
    setup_d = json.load(fjson)

dp = {}
detR_cm = 2.54
dist_cm = 220.
geom_fac = np.pi*detR_cm**2/(4.*np.pi*dist_cm**2)
det_eff = 0.2
setup_d['setup']['Start time'] = 0  # Take all events
setup_d['setup']['End time'] = -1 # Take all events
for nshot in (101, 102, 103, 104, 105, 107):
    f_ha = '/shares/departments/AUG/users/git/DPSD/acq/1/%d/HA_%d.dat' %(nshot, nshot)
    setup_d['io']['HA*.dat file'] = f_ha

    dp[nshot] = dpsd_run.DPSD(setup_d)

for nshot in (101, 102, 103, 104, 105, 107):
    dt = dp[nshot].time_cnt[-1] - dp[nshot].time_cnt[0]
    total1 = np.sum(dp[nshot].cnt['neut1'])
    total2 = np.sum(dp[nshot].cnt['neut2'])
    rate = total2/dt
    source_estimate = rate/geom_fac/det_eff
    print(nshot)
    print('Total neut: %12.4e' %total1)
    print('Total neut with pile-up: %12.4e' %total2) 
    print('Neut rate %12.4e [1/s]' %rate)
    print('Source emission estimate %12.4e' %source_estimate)

# Plots

    plt.figure(nshot, (10,6))
    plt.subplot(1, 1, 1)
    plt.figtext(0.6, 0.7, 'Average rate %12.4e' %rate, ha='center')
    plt.figtext(0.6, 0.6, 'Source estimate %12.4e' %source_estimate, ha='center')
    plt.title('Count rates %d' %nshot)
    plt.xlabel('Time [s]')
    plt.ylabel('Neutron count [1/s]')
    plt.plot(dp[nshot].time_cnt, dp[nshot].cnt['neut1'] , label='Neut')
    plt.plot(dp[nshot].time_cnt, dp[nshot].cnt['gamma1'], label='Gamma')
    plt.plot(dp[nshot].time_cnt, dp[nshot].cnt['neut2'] , label='Neut2')
    plt.legend()
    plt.savefig('neut%d.pdf' %nshot)
plt.show()

