from __future__ import print_function
import sys
import numpy as np
from astropy.io import fits
from matplotlib import pyplot as plt
import seaborn as sns
from astropy.table import Table
from scipy.special import erfinv
from stats_utils import (get_stats_for_var, get_fuzzstats_for_var, 
                         tsq_from_tables)

IQR_OVER_SIGMA = 2.0*np.sqrt(2.0)*erfinv(0.5)
print('IQR / sigma =', IQR_OVER_SIGMA)
BINNINGS = [1, 2, 4, 8, 16, 32, 64, 128, 256]

def plot_tsq_vs_bright(ax, ion, binnings, colors, robust=False):
    assert len(binnings) == len(colors)
    for nbin, color in zip(binnings, colors):
        stats = get_stats_for_var('T', ion, nbin=nbin)
        fuzzstats0 = get_fuzzstats_for_var('Delta_T', ion, nbin=nbin, ifuzz=0)
        fuzzstats1 = get_fuzzstats_for_var('Delta_T', ion, nbin=nbin, ifuzz=1)
        assert np.all(stats['X Center'] == fuzzstats0['X Center'])
        log_S = stats['X Center']
        npix = stats['Count']/(nbin**2)
        weight = stats['Count']*10**log_S
        weight /= weight.sum()
        m = (npix >= 5) & (weight > 0.5/len(stats))
        tsq = tsq_from_tables(stats, stats, robust)
        tsq_err0 = tsq_from_tables(fuzzstats0, stats, robust)
        tsq_err1 = tsq_from_tables(fuzzstats1, stats, robust)
        ax.plot(log_S[m], tsq[m], color=color, label='{0} x {0}'.format(nbin))
        ax.plot(log_S[m], tsq_err0[m], color=color, lw=0.5, label=None)
        ax.plot(log_S[m], tsq_err1[m], color=color, lw=0.5, label=None)

def plot_corrected_tsq_vs_bright(ax, ion, binnings, colors, robust=False):
    assert len(binnings) == len(colors)
    for nbin, color in zip(binnings, colors):
        stats = get_stats_for_var('T', ion, nbin=nbin)
        fuzzstats0 = get_fuzzstats_for_var('Delta_T', ion, nbin=nbin, ifuzz=0)
        fuzzstats1 = get_fuzzstats_for_var('Delta_T', ion, nbin=nbin, ifuzz=1)
        assert np.all(stats['X Center'] == fuzzstats0['X Center'])
        log_S = stats['X Center']
        npix = stats['Count']/(nbin**2)
        weight = stats['Count']*10**log_S
        weight /= weight.sum()
        m = (npix >= 5) & (weight > 0.5/len(stats))
        tsq = tsq_from_tables(stats, stats, robust)
        tsq_err0 = tsq_from_tables(fuzzstats0, stats, robust)
        tsq_err1 = tsq_from_tables(fuzzstats1, stats, robust)
        tsq -= 0.5*(tsq_err0 + tsq_err1)
        ax.plot(log_S[m], tsq[m], color=color, label='{0} x {0}'.format(nbin))

if __name__ == '__main__':
    binnings = BINNINGS[:8]
    #colors = sns.dark_palette('orange',  n_colors=len(binnings))
    #colors = sns.color_palette('Paired',  n_colors=len(binnings))
    colors = sns.color_palette('husl',  n_colors=len(binnings))
    #colors = sns.color_palette('Set2',  n_colors=len(binnings))
    #colors = ['k']*len(binnings)
    for ion in 'N_II', 'S_III':
        for variant in '', '-robust':
            fig, ax = plt.subplots(1, 1)
            plot_tsq_vs_bright(ax, ion, binnings, colors,
                               robust=variant=='-robust')
            leg = ax.legend(ncol=2, fontsize='x-small',
                            title='Pixel Binning', frameon=True, fancybox=True)
            leg.get_title().set_fontsize('x-small')
            ax.set_yscale('log')
            ax.set_ylim(1e-6, 3.0)
            ax.set_xlabel(r'$\log_{{10}}$( {} )'.format(WAVS[ion]))
            ax.set_ylabel('Partial $t^{{\,2}}$([{}])'.format(ion.replace('_', ' ')))
            fig.set_size_inches(4, 4)
            figname = 'muse-tsq{}-vs-bright-{}.pdf'.format(variant, ion)
            fig.tight_layout()
            fig.savefig(figname)
            print(figname)

            fig, ax = plt.subplots(1, 1)
            plot_corrected_tsq_vs_bright(ax, ion, binnings, colors,
                                         robust=variant=='-robust')
            leg = ax.legend(ncol=2, fontsize='x-small',
                            title='Pixel Binning', frameon=True, fancybox=True)
            leg.get_title().set_fontsize('x-small')
            ax.set_yscale('log')
            ax.set_ylim(1e-6, 3.0)
            ax.set_xlabel(r'$\log_{{10}}$( {} )'.format(WAVS[ion]))
            ax.set_ylabel(r'Noise-Corrected Partial $t^{{\,2}}$([{}])'.format(ion.replace('_', ' ')))
            fig.set_size_inches(4, 4)
            figname = 'muse-corrected-tsq{}-vs-bright-{}.pdf'.format(variant, ion)
            fig.tight_layout()
            fig.savefig(figname)
            print(figname)
