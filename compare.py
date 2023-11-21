import os, json
import numpy as np
import dpsd_run
import matplotlib.pylab as plt

f_json = 'settings/default.json'
with open(f_json) as fjson:
    setup_d = json.load(fjson)

nshot = 29795
setup_d['io']['Shots'] = nshot

#2nd harmonic

if nshot == 29783:
    t_nb = [67, 69]
    t_rf = [71, 73]
    t_shift = 65.85

if nshot == 29795:
    t_nb = [69, 71]
    t_rf = [72, 74]
    t_shift = 67.75

if nshot == 30669:
#        t_nb = [92, 92.8]
#        t_rf = [93, 93.8]
    t_nb = [95.9, 96.7]
    t_rf = [97.2, 97.6]
    t_shift = 90.85

if nshot == 30670: # very weak
    t_nb = [67, 67.7]
    t_rf = [68.05, 68.25]
    t_shift = 65.85

if nshot == 30672:
    t_shift = 65.85
    if phase == 1: # NBI 3
        t_nb = [ 67, 67.8]
        t_rf = [67.9, 68.8]
    if phase == 2: # NBI 3+8
        t_nb = [69.95, 70.6]
        t_rf = [ \
              [68.88, 68.97], \
              [69.11, 69.16], \
              [69.42, 69.48], \
              [69.58, 69.73] ]
    if phase == 3: # NBI 3+6
        t_nb = [70.8, 71.6]
        t_rf = [71.8, 72.6]

if nshot == 31589:
    t_shift = 65.35
    t_down1 = [ \
        [67.55, 67.58], \
        [67.72, 67.75], \
        [67.92, 67.96], \
        [68.14, 68.17] ]
    t_up1 = [ \
        [67.61, 67.7], \
        [67.8, 67.88], \
        [68, 68.1] ]

    t_rf = [68.5, 68.8]
    t_nb   = [69.5, 71.2]
#        t_nb = t_down1
#        t_rf = t_up1

if nshot == 31590:
    t_shift = 65.35
    if phase == 1: # NBI 3
        t_nb = [66.2, 67.2]
        t_rf = [67.4, 68.2]
    if phase == 2: # NBI 3 + 1
        t_nb = [ \
              [69.45, 69.75], \
              [69.9, 70.2] ]
        t_rf = [68.65, 69.15]
    if phase == 3: # NBI 3 + 4
        t_nb = [70.5, 71.2]
        t_rf = [ \
              [71.55, 71.8], \
              [71.9, 72.15] ]

#3rd

if nshot == 32573:
    t_nb = [2.85, 3.21]
    t_rf = [3.522, 4.22]

if nshot == 32604:
    t_nb = [4.86, 5.01]
#        t_rf = [4.32, 4.8]
    t_rf = [3.3, 3.48]

if nshot == 32605:
    t_nb = [2.46, 2.8]
    t_rf = [2.1, 2.4]

if nshot == 34227:
    t_nb = [8.05, 8.2]
    t_rf = [8.35, 8.65]

if nshot == 35448:
    t_shift = 0.05
    t_nb = [2.8, 3.8]
    t_rf = [4.2, 5.5]

if nshot == 35492:
    t_shift = 0.05
    t_nb = [2.5, 3.]
    t_rf = [4., 4.7]

if nshot == 35493:
    t_shift = 0.05
    t_nb = [2., 2.4]
    t_rf = [4.2, 4.6]

if nshot == 36557:
    t_shift = 0.05
#        t_rf = [3.14, 3.44]   # 2nd
    t_rf = [4.15, 4.45]   # 3rd
    t_nb = [3.6, 3.95]    # No RF

if nshot == 37949:
    t_shift = 0.05
    t_nb = [5.1, 5.4]
    t_rf = [5.6, 5.9]

# Compare coils on, coils off

if nshot in (34570, 34571):
    t_nb = [ \
         [1.67, 1.8], \
         [2.05, 2.20], \
         [2.45, 2.60], \
         [2.83, 3.0] ]
    t_rf = [ \
         [1.83, 1.99], \
         [2.24, 2.40], \
         [2.65, 2.80], \
         ]


nb = dpsd_run.DPSD(setup_d, t_ranges=t_nb)
rf = dpsd_run.DPSD(setup_d, t_ranges=t_rf)

plt.figure('DPSD compare', figsize=(12, 5))
plt.subplot(1, 2, 1)
plt.plot(nb.phs['neut1'], 'b-', label='NBI')
plt.plot(rf.phs['neut1'], 'r-', label='NBI&RF')
fac = np.mean(rf.phs['neut1'][50:150])/np.mean(nb.phs['neut1'][50:150])
plt.plot(fac*nb.phs['neut1'], 'b--', label='%5.2fxNBI' %fac)
plt.xlim([0, 400])
plt.ylim([0, 3000])
plt.legend()

plt.subplot(1, 2, 2)
plt.semilogy(nb.phs['neut1'], 'b-', label='NBI')
plt.semilogy(rf.phs['neut1'], 'r-', label='NBI&RF')
fac = np.mean(rf.phs['neut1'][50:150])/np.mean(nb.phs['neut1'][50:150])
plt.semilogy(fac*nb.phs['neut1'], 'b--', label='%5.2fxNBI' %fac)
plt.xlim([0, 400])
plt.ylim([1, 3000])
plt.legend()

plt.show()
