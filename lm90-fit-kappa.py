from __future__ import print_function
dataM=[["", "2.5e4", "3.2e4", "6.4e4", "1.28e5", "2.56e5", "5.12e5", "1.17e6", "<- T"], ["\\xi", "1e-3", "1.05e-3", "1.3e-3", "2.1e-3", "2.6e-3", "3.2e-3", "3.8e-4", "<- Kn"], [2.5, 0.99, 1.0, 1.0, 1.0, 1.0, 1.0, 0.99, ""], [3, 0.99, 0.99, 0.99, 1.0, 1.0, 1.0, 0.97, ""], [3.5, 0.99, 1.0, 1.01, 1.03, 1.06, 1.09, 0.93, ""], [4, 1.03, 1.04, 1.1, 1.24, 1.38, 1.59, 0.86, ""], [4.5, 1.2, 1.22, 1.49, 2.2, 3.03, 4.7, 0.78, ""], [5, 1.89, 1.96, 3.39, 9.01, 20.0, 39.0, 0.71, ""], [5.5, 4.96, 4.96, 25.6, "1.60e2", "5.56e2", "7.99e2", 0.68, ""], [6, 32.3, 43.7, "1.36e3", "1.21e4", "4.39e4", "2.84e4", 0.68, ""]]
dataB=[["", "2.5e4", "3.2e4", "6.4e4", "1.28e5", "2.56e5", "<- T"], ["\\xi", "2.2e-4", "6.2e-4", "5.9e-3", "2.8e-2", "7.9e-2", "<- Kn"], [2.5, 1.0, 1.01, 1.03, 1.16, 1.17, ""], [3, 1.0, 1.03, 1.23, 2.37, 2.3, ""], [3.5, 1.0, 1.12, 2.92, 12.4, 8.97, ""], [4, 1.01, 1.5, 32.3, "1.50e2", 67.7, ""], [4.5, 1.01, 12.7, "1.18e3", "3.31e3", "9.64e2", ""], [5, 1.02, "1.15e3", "1.02e5", "1.55e5", "2.66e4", ""], [5.5, 3.64, "1.37e5", "8.07e6", "6.38e6", "9.37e5", ""], [6, 85.3, "1.66e7", "3.34e9", "1.79e9", "1.40e8", ""]]
import sys
import numpy as np
from matplotlib import pyplot as plt
from astropy.table import Table
import seaborn as sns
from kappa_utils import f_M, f_CH, f_kappa

tab1 = Table(names=dataM[1], rows=dataM[2:])
tab2 = Table(names=dataB[1], rows=dataB[2:])


energy = np.logspace(-2, 2, 500)

sns.set_palette('hls', 7)
fig, ax = plt.subplots(1, 1)
ax.plot(energy, 1e7*f_M(energy),
        lw=7, alpha=0.1, color='k', label='Maxwellian, $10^{7} f_M$')
for kappa in 5.0, 10.0, 20.0, 40.0, 80.0, 160.0, 320.0:
    ax.plot(energy, f_kappa(energy, kappa)/f_M(energy),
            lw=1, alpha=0.5, label=r'$\kappa = {:.1f}$'.format(kappa))

for Kn in tab1.colnames[1:-2]:
    ax.plot(tab1[r'\xi']**2, tab1[Kn], '--', lw=3,
            label='McWhirter, Kn = {}'.format(Kn))

for Kn in tab2.colnames[1:-2]:
    ax.plot(tab2[r'\xi']**2, tab2[Kn], '-.', lw=3,
            label='Barlow, Kn = {}'.format(Kn))

ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlim(1.0, 100.0)
ax.set_ylim(0.1, 3e10)
ax.legend(fontsize='x-small', loc='best', ncol=2)
ax.set_xlabel(r'$E\, /\, k T$')
ax.set_ylabel(r'Excess over Maxwellian: $f\, /\, f_M$')
ax.set_title('Comparison of Ljepojevic & Burgess (1990) with kappa distributions')
figname = sys.argv[0].replace('.py', '.pdf')
fig.set_size_inches(7, 7)
fig.tight_layout()
fig.savefig(figname)
print(figname)
