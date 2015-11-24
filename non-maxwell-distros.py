from __future__ import print_function
import sys
import numpy as np
from matplotlib import pyplot as plt
import seaborn as sns
from kappa_utils import f_M, f_CH, f_kappa

energy = np.logspace(-2, 2, 500)

fig, ax = plt.subplots(1, 1)
ax.plot(energy, 1e7*f_M(energy), lw=7, alpha=0.1, color='k', label='Maxwellian, $10^{7} f_M$')
for kappa in 1.75, 3.0, 5.0, 10.0, 20.0, 100.0:
    ax.plot(energy, f_kappa(energy, kappa)/f_M(energy), lw=3, alpha=0.5, label=r'$\kappa = {:.1f}$'.format(kappa))
for a, b in (2, 0.05), (5, 0.05), (18, 0.2):
    ax.plot(energy, f_CH(energy, a, b)/f_M(energy), ls='--', lw=1.5, label='$T_C/T_H = {}$; $n_C/n_H = {:.2f}$'.format(int(a), b))

ax.set_xscale('log')
ax.set_yscale('log')
ax.set_ylim(0.1, 3e7)
ax.legend(fontsize='small', loc='middle left', ncol=2)
ax.set_xlabel(r'$E\, /\, k T$')
ax.set_ylabel(r'Excess over Maxwellian: $f\, /\, f_M$')
figname = sys.argv[0].replace('.py', '.pdf')
fig.set_size_inches(7, 5)
fig.tight_layout()
fig.savefig(figname)
print(figname)
