from __future__ import print_function
import sys
import numpy as np
from astropy.io import fits
from astropy.table import Table, Column
from matplotlib import pyplot as plt
import seaborn as sns
import pyregion
t2dir = '/Users/will/Work/RubinWFC3/Tsquared'
sys.path.append(t2dir)
from misc_utils import sanitize_string
from stats_utils import stats_vs_x

def prepare_string(label):
    newlabel = label.replace('log10 (N\' / N)', 'log10-n')
    return sanitize_string(newlabel)

try:
    fuzz = sys.argv[1]
    binning = sys.argv[2]
    region = sys.argv[3]
except IndexError:
    sys.exit('Usage: {} FUZZ BINNING REGION'.format(sys.argv[0]))

file_patterns = {
    '[N II] Delta T'  : 'muse-delta-Te-fuzz{}-bin{}.fits',
    '[S II] log10 (N\' / N)'  : 'muse-log10-n-fuzz{}-bin{}.fits',
    '[S III] Delta T' : 'muse-delta-Te-iii-fuzz{}-bin{}.fits',
    '[Cl III] log10 (N\' / N)': 'muse-log10-n-iii-fuzz{}-bin{}.fits',
    'S(Pa 9)'    : 'linesum-H_I-9229-fuzz{}-bin{}.fits',
    'S(6716)'    : 'linesum-S_II-6716-fuzz{}-bin{}.fits',
    'S(6583)'    : 'linesum-N_II-6583-fuzz{}-bin{}.fits',
    'S(9069)'    : 'linesum-S_III-9069-fuzz{}-bin{}.fits',
    'S(5518)'    : 'linesum-Cl_III-5518-fuzz{}-bin{}.fits',
    'T([N II])'  : 'muse-derived-Te-fuzz{}-bin{}.fits',
    'N([S II])'  : 'muse-derived-Ne-fuzz{}-bin{}.fits',
    'T([S III])' : 'muse-derived-Te-iii-fuzz{}-bin{}.fits',
    'N([Cl III])': 'muse-derived-Ne-iii-fuzz{}-bin{}.fits',
}
minmax = {
    '[N II] Delta T'  : [-5000.0, 5000.0], 
    '[S II] log10 (N\' / N)'  : [-1.5, 1.5], 
    '[S III] Delta T' : [-5000.0, 5000.0], 
    '[Cl III] log10 (N\' / N)': [-1.5, 1.5], 
    'S(Pa 9)'    : [4.0e3, 4.0e5],
    'S(5755)'    : [600.0, 6.0e5],
    'S(6716)'    : [600.0, 6.0e5],
    'S(6312)'    : [600.0, 6.0e5],
    'S(5518)'    : [600.0, 6.0e5],
    'T([N II])'  : [5000.0, 15000.],
    'N([S II])'  : [100.0, 1e5],
    'T([S III])' : [5000.0, 15000.], 
    'N([Cl III])': [100.0, 1e5],
}



hdus = {k: fits.open('Linemaps/' + v.format(fuzz, binning))['SCALED']
        for k, v in file_patterns.items()}

pairs = [
    ['S(Pa 9)', '[N II] Delta T'],
    ['S(Pa 9)', '[S II] log10 (N\' / N)'],
    ['S(Pa 9)', '[S III] Delta T'],
    ['S(Pa 9)', '[Cl III] log10 (N\' / N)'],
]

# Pairs that need to know about each other because of co-dependent ratios
complements = {
    '[S II] log10 (N\' / N)'   : 'T([N II])', 
    '[N II] Delta T'   : 'N([S II])', 
    '[Cl III] log10 (N\' / N)' : 'T([S III])', 
    '[S III] Delta T'  : 'N([Cl III])', 
}

weighting_maps = {
    '[N II] Delta T'  : 'S(6583)', 
    '[S II] log10 (N\' / N)'  : 'S(6716)', 
    '[S III] Delta T' : 'S(9069)', 
    '[Cl III] log10 (N\' / N)': 'S(5518)', 
}

mm = np.isfinite(hdus['S(Pa 9)'].data) & (hdus['S(Pa 9)'].data > 0.0) 
if 'sweet' in region.lower():
    title = 'Orion S'
    include = pyregion.open(t2dir + '/will-nii-sweet-spot-wcs.reg')
    exclude = pyregion.open(t2dir + '/will-nii-exclude-new.reg')
    mm = mm & include.get_mask(hdu=hdus['S(Pa 9)']) & (~exclude.get_mask(hdu=hdus['S(Pa 9)']))

    include = pyregion.open(t2dir + '/will-sii-sweet-spot-wcs.reg')
    exclude = pyregion.open(t2dir + '/will-sii-exclude-new.reg')
    mm = mm & include.get_mask(hdu=hdus['S(Pa 9)']) & (~exclude.get_mask(hdu=hdus['S(Pa 9)']))
else:
    title = 'Full nebula'

for xlabel, ylabel in pairs:
    hduA = hdus[xlabel]
    hduB = hdus[ylabel]
    xmin, xmax = minmax[xlabel]
    ymin, ymax = minmax[ylabel]

    m = mm & np.isfinite(hduA.data) & np.isfinite(hduB.data)
    m = m & (hduA.data > xmin) & (hduB.data > ymin) 
    m = m & (hduA.data < xmax) & (hduB.data < ymax) 
    m = m & (hduA.data != 0.0) & (hduB.data != 0.0) 

    wlabels = []
    for label in xlabel, ylabel:
        # Make a list of quantities to use as weights for this map
        if label in weighting_maps:
            wlabels.append(weighting_maps[label])
        else:
            wlabels.append(label)
        # Also check complementary vars for being out-of-range
        if label in complements:
            zlabel = complements[label]
            zmin, zmax = minmax[zlabel]
            hduC = hdus[zlabel]
            m = m & (hduC.data > zmin) & (hduC.data < zmax)

    x = hduA.data[m]
    y = hduB.data[m]
    w = np.zeros_like(x)
    for wlabel in set(wlabels):
        w += hdus[wlabel].data[m]

    gamma = 1.5
    # Red tinted colormap for the error histograms
    cmap = sns.light_palette((30, 50, 30), input="husl", as_cmap=True)
    fig, ax = plt.subplots(1, 1)
    # Save a sorted, linear version of the weights
    ww = np.sort(w)
    # Use a log scale for surface brightnesses and densities
    if xlabel.startswith('S') or xlabel.startswith('N'):
        x = np.log10(x)
        xmin, xmax = np.log10(xmin), np.log10(xmax)
        xlabel = 'log10[ {} ]'.format(xlabel)
    if ylabel.startswith('S') or ylabel.startswith('N'):
        y = np.log10(y)
        ymin, ymax = np.log10(ymin), np.log10(ymax)
        ylabel = 'log10[ {} ]'.format(ylabel)
    H, xedges, yedges = np.histogram2d(x, y, 51,
                                       [[xmin, xmax], [ymin, ymax]],
                                       weights=w)
    stats_table = stats_vs_x(x, y, xedges)
    ax.imshow((H.T)**(1.0/gamma), extent=[xmin, xmax, ymin, ymax],
              interpolation='nearest', aspect='auto', origin='lower', 
              cmap=cmap, alpha=1.0)

    ax.set_xlim(xmin, xmax)
    ax.set_ylim(ymin, ymax)
    ax.set_xlabel(xlabel, fontsize='large')
    ax.set_ylabel(ylabel, fontsize='large')
    ax.grid(color='w', lw=0.5, alpha=0.5, zorder=100)

    fig.set_size_inches(3, 3)
    fig.tight_layout()
    pattern = 'muse-{}-{}-histogram-{}-fuzz{}-bin{}.pdf'
    pltfile = pattern.format(prepare_string(xlabel), prepare_string(ylabel),
                             region, fuzz, binning)

    fig.savefig(pltfile)
    statsfile = pltfile.replace('.pdf', '.tab').replace('-histogram-', '-stats-')
    stats_table.write(statsfile, format='ascii.tab')

    histfile = statsfile.replace('-stats-', '-hist-1d-')

    # Now save the full histogram for each tertile in the CDF of weights
    wcdf = np.cumsum(ww)/np.sum(ww)
    w1, w2 = [np.max(ww[100*wcdf < percentile]) for percentile in [33, 67]]
    print('Tertiles of CDF of weights:')
    print(w1, w2)
    m1 = (w <= w1)
    m2 = (w > w1) & (w <= w2)
    m3 = (w > w2) 
    HH1, yedges = np.histogram(y[m1], bins=51, range=[ymin, ymax], weights=w[m1])
    HH2, yedges = np.histogram(y[m2], bins=51, range=[ymin, ymax], weights=w[m2])
    HH3, yedges = np.histogram(y[m3], bins=51, range=[ymin, ymax], weights=w[m3])

    ycenters = 0.5*(yedges[:-1] + yedges[1:])

    hist_table = Table()
    hist_table.add_columns([
        Column(name=ylabel, data=ycenters),
        Column(name='T1', data=HH1),
        Column(name='T2', data=HH2),
        Column(name='T3', data=HH3),
    ])
    hist_table.write(histfile, format='ascii.tab')

    print(pltfile)
