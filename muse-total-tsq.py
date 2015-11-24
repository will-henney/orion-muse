from __future__ import print_function
import sys
import numpy as np
from astropy.io import fits
from astropy.table import Table, join
from matplotlib import pyplot as plt
import seaborn as sns
from stats_utils import (get_stats_for_var,
                         get_fuzzstats_for_var,
                         tsq_from_tables)

BINNINGS = [1, 2, 4, 8, 16, 32, 64, 128, 256]
 
def nanaverage(a, w):
    '''Weighted average of `a` with weights `w`.  Just like np.average,
but ignoring NaNs'''
    return np.nansum(a*w)/np.nansum(w)
    

def total_combined_tsq(stats, robust=False, meanstats=None):
    '''Find weighted tsq over all brightnesses, including contribution
from trend in mean T with brightness.  `stats` is a table of
temperature statistics as a function of brightness bin.  If `robust` is
`True` then median and IQR are used instead of mean and std.  If
`stats` gives statistics for delta T instead of T (for instance, from
fuzzing), then 'meanstats' must be supplied, giving the equivalent
stats for T.

    '''
    if meanstats is None:
        tsq_vs_bright = tsq_from_tables(stats, stats, robust)
    else:
        checkbins = np.all(stats['X Center'] == meanstats['X Center'])
        assert checkbins, 'Error: brightness bins do not coincide'
        tsq_vs_bright = tsq_from_tables(stats, meanstats, robust)
    # Weight each brightness bin by the total flux from that bin
    weights = stats['Count']*(10**stats['X Mean'])
    average_tsq = nanaverage(tsq_vs_bright, weights)
    # Contribution from trend of mean (or median) T with brightness
    av_type = 'Y Median' if robust else 'Y Mean'
    T0_vs_bright = stats[av_type]
    average_T0 = nanaverage(T0_vs_bright, weights)
    average_T00 = average_T0
    if meanstats is not None:
        # Case where T0 is actually a Delta T, so we need to add on
        # the "real" mean T
        average_T00 += nanaverage(meanstats[av_type], weights)
    trend_tsq_vs_bright = ((T0_vs_bright - average_T0)/average_T00)**2
    
    trend_tsq = nanaverage(trend_tsq_vs_bright, weights)
    return {
        't^2': average_tsq + trend_tsq,
        'T': average_T0,
        'trend fraction': trend_tsq / (average_tsq + trend_tsq),
    }


def find_tsq_vs_scale(ion, robust=False, ifuzz=0, region='full'):
    '''Calculate total t^2 versus binning. Returns an astropy.table'''
    results = {'nbin': BINNINGS, 't^2': [], 'T': [], 'trend fraction': []}
    fuzz_results = {'nbin': BINNINGS, 't^2': [], 'T': [], 'trend fraction': []}
    for nbin in BINNINGS:
        stats = get_stats_for_var('T', ion,
                                  region=region, nbin=nbin)
        fstats = get_fuzzstats_for_var('Delta_T', ion,
                                       region=region, nbin=nbin, ifuzz=ifuzz)
        for k, v in total_combined_tsq(stats, robust).items():
            results[k].append(v)
        for k, v in total_combined_tsq(fstats, robust, meanstats=stats).items():
            fuzz_results[k].append(v)
    tab = Table(results)
    ftab = Table(fuzz_results)
    for t in tab, ftab:
        t['T'].format = '{:.0f}'
        t['t^2'].format = '{:.5f}'
        t['trend fraction'].format = '{:.2f}'
    return join(tab, ftab, keys='nbin', table_names=['', 'd '],
                uniq_col_name='{table_name}{col_name}')



for ion in 'N_II', 'S_III':
    for variant in '', '-robust':
        for region in 'full', 'sweet':
            tab = find_tsq_vs_scale(ion, robust=(variant=='-robust'), region=region) 
            fn = 'muse-total-tsq-{}{}-{}.tab'.format(ion, variant, region)
            tab.write(fn, format='ascii.tab')
            print(fn)
#            print('[[{}]]'.format(fn))
