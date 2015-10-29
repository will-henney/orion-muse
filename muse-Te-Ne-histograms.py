from __future__ import print_function
import sys
import numpy as np
from astropy.io import fits
from matplotlib import pyplot as plt
import pyregion
t2dir = '/Users/will/Work/RubinWFC3/Tsquared'
sys.path.append(t2dir)
from pyneb_utils import rsii_T_den, rnii_T_den
from misc_utils import sanitize_string


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
    'S(6716)'    : 'linesum-S_II-6716-bin{}.fits',
    'S(6312)'    : 'linesum-S_III-6312-bin{}.fits',
    'S(5518)'    : 'linesum-Cl_III-5518-bin{}.fits',
}

minmax = {
    'T([N II])'  : [5000.0, 15000.],
    'N([S II])'  : [100.0, 1e5],
    'T([S III])' : [5000.0, 15000.], 
    'N([Cl III])': [100.0, 1e5],
    'S(Pa 9)'    : [0.01, 8.0],
    'S(5755)'    : [0.01, 8.0],
    'S(6716)'    : [0.01, 8.0],
    'S(6312)'    : [0.01, 8.0],
    'S(5518)'    : [0.01, 8.0],
}



hdus = {k: fits.open('Linemaps/' + v.format(suffix))['SCALED']
        for k, v in file_patterns.items()}

# Rescale brightnesses to give nicer numbers
for k in hdus:
    if k.startswith('S'):
        hdus[k].data /= 1.e5

pairs = [
    ['T([N II])', 'T([S III])'],
    ['N([S II])', 'N([Cl III])'],
    ['N([S II])', 'T([N II])'],
    ['N([Cl III])', 'T([S III])'],
    # ['S(Pa 9)', 'N([S II])'],
    # ['S(Pa 9)', 'N([Cl III])'],
    # ['S(Pa 9)', 'T([N II])'],
    # ['S(Pa 9)', 'T([S III])'],
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
    'T([N II])': 'S(5755)', 
    'N([Cl III])': 'S(5518)', 
    'T([S III])': 'S(6312)', 
}

mm = np.isfinite(hdus['S(Pa 9)'].data) & (hdus['S(Pa 9)'].data > 0.0) 
if 'sweet' in region.lower():
    title = 'Orion S'
    include = pyregion.open(t2dir + '/will-nii-sweet-spot-wcs.reg')
    exclude = pyregion.open(t2dir + '/will-nii-exclude-wcs.reg')
    mm = mm & include.get_mask(hdu=hdus['S(Pa 9)']) & (~exclude.get_mask(hdu=hdus['S(Pa 9)']))

    include = pyregion.open(t2dir + '/will-sii-sweet-spot-wcs.reg')
    exclude = pyregion.open(t2dir + '/will-sii-exclude-wcs.reg')
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
            m = m & (hduC.data > zmin) & (hduC.data < zmax)

    x = hduA.data[m]
    y = hduB.data[m]
    w = np.zeros_like(x)
    for wlabel in set(wlabels):
        w += hdus[wlabel].data[m]
    gamma = 1.5
    fig, ax = plt.subplots(1, 1)
    # Use a log scale for surface brightnesses and densities
    if xlabel.startswith('S') or xlabel.startswith('N'):
        x = np.log10(x)
        xmin, xmax = np.log10(xmin), np.log10(xmax)
        xlabel = 'log10[ {} ]'.format(xlabel)
    if ylabel.startswith('S') or ylabel.startswith('N'):
        y = np.log10(y)
        ymin, ymax = np.log10(ymin), np.log10(ymax)
        ylabel = 'log10[ {} ]'.format(ylabel)
    H, xedges, yedges = np.histogram2d(x, y, 50,
                                       [[xmin, xmax], [ymin, ymax]],
                                       weights=w)
    ax.imshow((H.T)**(1.0/gamma), extent=[xmin, xmax, ymin, ymax],
              interpolation='nearest', aspect='auto', origin='lower', 
              cmap=plt.cm.gray_r, alpha=1.0)

    ax.set_xlim(xmin, xmax)
    ax.set_ylim(ymin, ymax)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.grid()

    fig.set_size_inches(4, 4)
    fig.tight_layout()
    pltfile = 'muse-{}-{}-histogram-{}{}.pdf'.format(sanitize_string(xlabel),
                                                     sanitize_string(ylabel),
                                                     region, suffix)

    fig.savefig(pltfile)
    print(pltfile)
