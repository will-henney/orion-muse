from __future__ import print_function
import sys
import numpy as np
from matplotlib import pyplot as plt
import seaborn as sns
from astropy.table import Table

TERTILES = ['T1', 'T2', 'T3']
BINNINGS = [1, 4, 16, 64]

LINETYPES = ['-', '--', '-.', ':']
LINEWIDTHS = [1, 1.5, 2.0, 2.5]
# Sweet spot doesn't have enough pixels to tolerate 64x64 binning
MAXBINNING = {'full': 64, 'sweet': 16}

quant_label = {
    'N_II_Delta_T'  : '[N II] Delta T',
    'S_III_Delta_T' : '[S III] Delta T',
    'S_II_log10-n'  : '[S II] log10 (N\' / N)',
    'Cl_III_log10-n': '[Cl III] log10 (N\' / N)',
}

quant_weight_text = {
    'N_II_Delta_T': '[N II] flux',
    'S_III_Delta_T': '[S III] flux',
    'S_II_log10-n': '[S II] flux',
    'Cl_III_log10-n': '[Cl III] flux',
}

def get_data(xlabel='log10_S_Pa_9_', ylabel='T_N_II', binning=1, region='full', ifuzz=0):
    tabfile = 'muse-{xlabel}-{ylabel}-hist-1d-{region}-fuzz{ifuzz:03d}-bin{binning:03d}.tab'.format(
        xlabel=xlabel, ylabel=ylabel, binning=binning, region=region, ifuzz=ifuzz)
    return Table.read(tabfile, format='ascii.tab')


if __name__ == '__main__':
    centiles = TERTILES
    try: 
        quant = sys.argv[1]
        region = sys.argv[2]
  
    except IndexError:
        print('Usage: {} QUANTITY ("full"|"sweet")'.format(sys.argv[0]))

    clabels = 'Brightest pixels', 'Mid-brightness pixels', 'Faintest pixels'
    sns.set_palette('dark')
    fig, axes = plt.subplots(len(centiles), 1, sharex=True, sharey=True)
    for binning, ls, lw in zip(BINNINGS, LINETYPES, LINEWIDTHS):
        if binning > MAXBINNING[region]:
            continue
        for ifuzz in range(10):
            data = get_data(binning=binning, region=region, ylabel=quant, ifuzz=ifuzz)
            if ifuzz == 0:
                sumdata = data
            else:
                for centile in centiles:
                    sumdata[centile] += data[centile]
        x = sumdata[quant_label[quant]]

        for ax, centile in zip(axes, centiles[::-1]):
            y = sumdata[centile]
            # x = 0.5*(x[:-1:2] + x[1::2])
            # y = 0.5*(y[:-1:2] + y[1::2])
            y /= y.sum()
            ax.plot(x, y, lw=lw, ls=ls,
                    label='{0} x {0}'.format(binning))
            ax.fill_between(x, y, alpha=0.2, color=(0.7, 0.3, 0.3))
    for ax, clabel in zip(axes, clabels):
        ax.set_yscale('log')
        ax.set_ylim(1e-3, 0.5)
        ax.text(0.02, 0.95, clabel, transform=ax.transAxes, fontsize='small')
    axes[-1].set_xlabel(quant_label[quant])
    axes[1].set_ylabel('Fraction of ' + quant_weight_text[quant])
    legend = axes[-1].legend(fontsize='small', title='Spatial binning')
    legend.get_title().set_fontsize('small')
    fn = 'muse-plot-{}-noise-hist1d-{}.pdf'.format(quant, region)
    fig.set_size_inches(6, 9)
    fig.tight_layout()
    fig.savefig(fn)
