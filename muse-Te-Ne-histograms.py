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


try:
    region = sys.argv[1]
    suffix = sys.argv[2]
except IndexError:
    sys.exit('Usage: {} REGION [SUFFIX]'.format(sys.argv[0]))

file_patterns = {
    'T([N II])'  : 'muse-derived-Te-bin{}.fits',
    'N([S II])'  : 'muse-derived-Ne-bin{}.fits',
    'T([S III])' : 'muse-derived-Te-iii-bin{}.fits',
    'N([Cl III])': 'muse-derived-Ne-iii-bin{}.fits',
    'S(Pa 9)'    : 'linesum-H_I-9229-bin{}.fits',
    'S(5755)'    : 'linesum-N_II-5755-bin{}.fits',
    'S(6583)'    : 'linesum-N_II-6583-bin{}.fits',
    'S(6716)'    : 'linesum-S_II-6716-bin{}.fits',
    'S(6312)'    : 'linesum-S_III-6312-bin{}.fits',
    'S(5518)'    : 'linesum-Cl_III-5518-bin{}.fits',
}

minmax = {
    'T([N II])'  : [5000.0, 15000.],
    'N([S II])'  : [100.0, 1e5],
    'T([S III])' : [5000.0, 15000.], 
    'N([Cl III])': [100.0, 1e5],
    'S(Pa 9)'    : [4.0e3, 4.0e5],
    'S(5755)'    : [600.0, 6.0e5],
    'S(6716)'    : [600.0, 6.0e5],
    'S(6312)'    : [600.0, 6.0e5],
    'S(5518)'    : [600.0, 6.0e5],
}



hdus = {k: fits.open('Linemaps/' + v.format(suffix))['SCALED']
        for k, v in file_patterns.items()}

true_shape = hdus['S(Pa 9)'].data.shape
for lab, hdu in hdus.items():
    assert hdu.data.shape == true_shape, '{} has shape {}, should be {}'.format(hdu.fileinfo()['file'], hdu.data.shape, true_shape)

pairs = [
    ['T([N II])', 'T([S III])'],
    ['N([S II])', 'N([Cl III])'],
    ['N([S II])', 'T([N II])'],
    ['N([Cl III])', 'T([S III])'],
    ['S(Pa 9)', 'N([S II])'],
    ['S(Pa 9)', 'N([Cl III])'],
    ['S(Pa 9)', 'T([N II])'],
    ['S(Pa 9)', 'T([S III])'],
    ['S(6716)', 'N([S II])'],
    ['S(5518)', 'N([Cl III])'],
    ['S(5755)', 'T([N II])'],
    ['S(6312)', 'T([S III])'],
]

# Pairs that need to know about each other because of co-dependent ratios
complements = {
    'N([S II])': 'T([N II])', 
    'T([N II])': 'N([S II])', 
    'N([Cl III])': 'T([S III])', 
    'T([S III])': 'N([Cl III])', 
}

weighting_maps = {
    'N([S II])': 'S(6716)', 
    'T([N II])': 'S(6583)', 
    'N([Cl III])': 'S(5518)', 
    'T([S III])': 'S(9069)', 
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
            print('Labels:', xlabel, ylabel, zlabel)
            print('Shapes:', hduA.data.shape, hduB.data.shape, hduC.data.shape)
            m = m & (hduC.data > zmin) & (hduC.data < zmax)

    x = hduA.data[m]
    y = hduB.data[m]
    w = np.zeros_like(x)
    for wlabel in set(wlabels):
        w += hdus[wlabel].data[m]
    gamma = 1.5
    # Blue tinted colormap for the data histograms
    cmap = sns.light_palette((260, 50, 30), input="husl", as_cmap=True)
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
    pltfile = 'muse-{}-{}-histogram-{}{}.pdf'.format(sanitize_string(xlabel),
                                                     sanitize_string(ylabel),
                                                     region, suffix)

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
