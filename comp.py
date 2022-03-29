import numpy as np
import dpsd, dixm
import matplotlib.pylab as plt

setup_d = dixm.DIX().xml2dict('/afs/ipp/home/g/git/python/neutrons/dpsd/xml/default.xml')
setup_d['Shot'] = 29795

t_nb = [22.2, 24.2]
t_rf = [24.77, 27.91]
t_shift = 67.75

nb = dpsd.DPSD(setup_d, t_ranges=t_nb)
rf = dpsd.DPSD(setup_d, t_ranges=t_rf)

plt.figure('DPSD compare', figsize=(12, 5))
plt.subplot(1, 2, 1)
plt.plot(nb.phs['neut1'], 'b-', label='NBI')
plt.plot(rf.phs['neut1'], 'r-', label='NBI&RF')
fac = np.mean(rf.phs['neut1'][50:150])/np.mean(nb.phs['neut1'][50:150])
plt.plot(fac*nb.phs['neut1'], 'b--', label='%5.2fxNBI' %fac)
plt.xlim([0, 400])
plt.ylim([0, 2500])
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
