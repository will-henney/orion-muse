import matplotlib.pyplot as plt
import pysynphot

filters = 'f547m', 'f469n', 'f673n'
MJDs = 'mjd#55233', 'mjd#55933', 'qyc'
linestyles = '-', '--', '-.'
for mjd, ls in zip(MJDs, linestyles):
    for f in filters:
        bp = pysynphot.ObsBandpass(','.join(['wfc3,uvis1', f, mjd]))
        plt.plot(bp.wave, bp.throughput, ls=ls, label='{} {}'.format(f, mjd))
plt.xlim(4500, 8000)
plt.legend(fontsize='small')
plt.savefig('wfc3-throughput-evolution.pdf')
