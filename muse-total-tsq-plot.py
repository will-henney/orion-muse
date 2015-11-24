from __future__ import print_function
import sys
import numpy as np
from matplotlib import pyplot as plt
import seaborn as sns
from astropy.table import Table

PIXEL_SCALE = 0.2
region_label = {'full': 'Full MUSE field', 'sweet': 'WFC3 field only'}


def line_plus_symbol_plot(ax, x, y, label=None, symbol='o'):
    '''Like normal ax.plot, but use both a solid line and symbols, all
with the same color

    '''
    color = next(ax._get_lines.color_cycle)
    ax.plot(x, y, color=color, label=label)
    ax.plot(x, y, symbol, color=color, label=None)


def plot_tsq_vs_scale(ion, region='full', variant='-robust', vscale=None):
    fn = 'muse-total-tsq-{}{}-{}.tab'.format(ion, variant, region)
    tab = Table.read(fn, format='ascii.tab')
    fig, ax = plt.subplots(1, 1)
    scale = tab['nbin']*PIXEL_SCALE
    line_plus_symbol_plot(ax, 2*scale, tab['t^2'], 'raw')
    line_plus_symbol_plot(ax, 2*scale, tab['d t^2'], 'noise')
    line_plus_symbol_plot(ax, 2*scale, tab['t^2'] - tab['d t^2'], 'corrected')
    if vscale is not None:
        ax.set_ylim(-0.02*vscale, 1.02*vscale)
    ax.set_xscale('log')
    ax.set_xlabel('Nyquist angular scale ($2 \\times$ bin width), arcsec')
    ax.set_ylabel('Total plane-of-sky $t^{{\,2}}$([{}])'
                  .format(ion.replace('_', ' ')))
    leg = ax.legend(fontsize='small', title=region_label[region],
                    frameon=True, fancybox=True)
    leg.get_title().set_fontsize('x-small')
    fig.set_size_inches(4, 4)
    figname = fn.replace('.tab', '.pdf')
    fig.tight_layout()
    fig.savefig(figname)
    print(figname)



sns.set_palette('husl',  n_colors=3)
for ion in 'N_II', 'S_III':
    for region in 'full', 'sweet':
        for variant in '-robust', '':
            vscale = 0.02 if ion == 'N_II' else 0.005
            plot_tsq_vs_scale(ion, region, variant, vscale)
