from __future__ import print_function
data_B=[["Kn", "kappa"], ["2.2e-4", 160], ["6.2e-4", 50], ["5.9e-3", 20], ["2.8e-2", 10]]
data_A=[["Kn", "kappa"], ["1e-3", 160], ["1.05e-3", 160], ["1.3e-3", 110], ["2.1e-3", 90], ["2.6e-3", 80], ["3.2e-3", 70]]
data_3=[["Region", "T / K", "n / pcc", "H", "ln \\Lambda", "\\lambda_{e}", "Kn", "\\kappa"], ["average loop", "3.2e6", "1.8e9", "3.33e-3", 21.26, "7.02e7", 0.3, 2], ["y=300-309", "3.2e6", "1.6e9", "3.33e-3", 21.32, "7.88e7", 0.34, 2]]
data_2=[["Region", "T / K", "n / pcc", "H", "ln \\Lambda", "\\lambda_{e}", "Kn", "\\kappa"], ["Coronal Hole", "2.5e4", "1.4e10", "6e3", 12.96, "9.04e2", 0.15, 13], ["Quiet Sun", "3.5e4", "1.8e9", "1.5e4", 14.49, "1.23e4", 0.82, 10], ["Active Region", "1e4", "1.3e10", "5.9e2", 11.62, "1.74e2", 0.29, 7]]
data_1=[["R/Rsun", "T / K", "n / pcc", "H", "ln \\Lambda", "\\lambda_{e}", "Kn", "n_{h}/n_{c}", "T_{h}/T_{c}", "\\kappa"], [1.0, "5e5", "3.8e8", 0.07, 19.26, "8.97e6", "1.8e-3", 0.05, "<2", 20], [1.25, "9e5", "1e7", 0.07, 21.96, "9.68e8", 0.2, 0.05, 5, 3], [1.5, "9e5", "1e6", 0.2, 23.11, "9.20e9", 0.66, "", "", -1], [2.0, "7e5", "2e5", 0.4, 23.54, "2.73e10", 0.98, "", "", -1], [2.4, "6e5", "1e5", 0.4, 23.65, "4.00e10", 1.44, 0.2, 18, 2]]
import sys
import numpy as np
from astropy.table import Table
from matplotlib import pyplot as plt
import seaborn as sns

def clean_data(data):
    """Replace empty strings with nans"""
    for row in data:
        for x in row:
            x = x or -1.0
    return data

d = {}
for label, data in (
        ['SolarWind-EE00', data_1],
        ['TR-DK11', data_2],
        ['CoronalLoop-D15', data_3],
        ['LB90-McWhirter', data_A],
        ['LB90-Barlow', data_B],
):
    cdata = clean_data(data)
    d[label] = Table(names=cdata[0], rows=cdata[1:])

sns.set_palette('hls', 7)
fig, ax = plt.subplots(1, 1)
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlim(1e-10, 10.0)
ax.set_ylim(1.0, 1.e6)
ax.set_xlabel(r'Knudsen number: $\mathsf{Kn} = \lambda/L$', fontsize='x-large')
ax.set_ylabel(r'Velocity distribution parameter: \large $\kappa$', fontsize='x-large')
ax.xaxis.set_ticks([1e-9, 1e-6, 1e-3, 1.0])

Kn = np.array([1e-12, 1.0, 100.0])
kappa1 = 1.5/np.sqrt(Kn)
kappa1[-1] = 1.5
kappa2 = 10/np.sqrt(Kn)
plt.fill_between(Kn, kappa1, kappa2, alpha=0.1, lw=0.0, color='k')

for dataset in d:
    Kn = d[dataset]['Kn'].astype('float')
    try:
        kappa = d[dataset]['kappa'].astype('float')
    except KeyError:
        kappa = d[dataset][r'\kappa'].astype('float')
    m = kappa > 0.0

    if 'LB90' in dataset:
        plotstyles = {'ls': '-', 'lw': 3}
    else:
        plotstyles = {'marker': 'o', 'ls': ''}
    plt.plot(Kn[m], kappa[m], label=dataset, **plotstyles)


plt.legend()
figfile = sys.argv[0].replace('.py', '.pdf')
fig.set_size_inches(6, 6)
fig.tight_layout()
fig.savefig(figfile)
print(figfile)
