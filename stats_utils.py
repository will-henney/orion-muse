import numpy as np
import scipy.stats as ss
from astropy.table import Table
from scipy.special import erfinv

WAVS = {'N_II': '5755', 'S_III': '6312',}
IQR_OVER_SIGMA = 2.0*np.sqrt(2.0)*erfinv(0.5)

def stats_vs_x(x, y, xedges):
    '''Calculate statistics column by column to complement a histogram

    Returns an astropy.table object
    '''
    rslt = {'Y Mean': [], 'Y Median': [], 'Y Std': [],
            'Y Lower Quartile': [], 'Y Upper Quartile': [],
            'Count': [], 'X Mean': []}
    for x1, x2 in zip(xedges[:-1], xedges[1:]):
        thisbin = (x >= x1) & (x <= x2)
        rslt['Y Mean'].append(y[thisbin].mean())
        rslt['Y Median'].append(np.median(y[thisbin]))
        rslt['Y Std'].append(y[thisbin].std())
        rslt['Count'].append(thisbin.sum())
        rslt['X Mean'].append(x[thisbin].mean())
        if thisbin.sum() > 0:
            rslt['Y Lower Quartile'].append(
                ss.scoreatpercentile(y[thisbin].ravel(), 25))
            rslt['Y Upper Quartile'].append(
                ss.scoreatpercentile(y[thisbin].ravel(), 75))
        else:
            rslt['Y Lower Quartile'].append(np.nan)
            rslt['Y Upper Quartile'].append(np.nan)
    for s in rslt:
        rslt[s] = np.array(rslt[s])
    rslt['X Center'] = 0.5*(xedges[:-1] + xedges[1:])
    return Table(rslt)


def get_stats_for_var(var='T', ion='N_II', region='full', nbin=1):
    xstring = 'log10_S_{}_'.format(WAVS[ion])
    ystring = '{}_{}'.format(var, ion)
    tabname = 'muse-{}-{}-stats-{}{:03d}.tab'.format(xstring, ystring, region, nbin)
    return Table.read(tabname, format='ascii.tab')

def get_fuzzstats_for_var(var='Delta_T', ion='N_II', region='full', nbin=1, ifuzz=0):
    xstring = 'log10_S_{}_'.format(WAVS[ion])
    ystring = '{}_{}'.format(ion, var)
    tabname = 'muse-{}-{}-stats-{}-fuzz{:03d}-bin{:03d}.tab'.format(
        xstring, ystring, region, ifuzz, nbin)
    return Table.read(tabname, format='ascii.tab')


def tsq_from_tables(tabvar, tabmean, robust=False):
    if robust:
        # Use Interquartile range as proxy for standard deviation
        sigma = (tabvar['Y Upper Quartile'].data
                 - tabvar['Y Lower Quartile'].data)/IQR_OVER_SIGMA
        # Use median as proxy for mean
        av = tabmean['Y Median'].data
    else:
        sigma = tabvar['Y Std'].data
        av = tabmean['Y Mean'].data

    return sigma**2/av**2
