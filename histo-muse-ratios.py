from __future__ import print_function
import sys
import numpy as np
from astropy.io import fits
from astropy.convolution import convolve, Gaussian2DKernel
from matplotlib import pyplot as plt
import seaborn as sns
import pyregion
from specplot1d_utils import plot_1d_spec_from_fits

plotpars = {
    '6716-6731': {'line1': '[S II] 6716', 'line2': '[S II] 6731',
                  'min': 0.4, 'max': 0.8},
    '6716-6731-N': {'line1': '[S II] 6716', 'line2': '[S II] 6731',
                  'min': 0.4, 'max': 0.8},
    '5755-6583': {'line1': '[N II] 5755', 'line2': '[N II] 6583',
                  'min': 0.00, 'max': 0.05},
    '4861-6563': {'line1': 'H I 4861', 'line2': 'H I 6563',
                  'min': 0.2, 'max': 0.35},
  }

titles_from_extra = {
    '': 'Fully continuum-corrected',
    '-naive': 'Uncorrected for continuum',
    '-flat': 'Continuum-corrected without color terms',
}

GAMMA = 1.0

cmap = sns.light_palette((260, 50, 30), input="husl", as_cmap=True)
# cmap = plt.cm.gray_r

def histogram_ratio_images(ratio_name, pars, extra='', regionfile=None):
    ratio_name_true = '-'.join(ratio_name.split('-')[:2])
    fn_true = 'LineMaps/ratio-{}.fits'.format(ratio_name_true)
    fn_syn = 'NebulioMUSE/synthetic{}-ratio-{}.fits'.format(extra, ratio_name)
    pltname = 'NebulioMUSE/synthetic{}-vs-true-calib-{}.pdf'.format(extra, ratio_name)
    hdu_true = fits.open(fn_true)[0]
    hdu_syn = fits.open(fn_syn)[0]
    # Flux of a strong line to weight the pixels
    hduf = fits.open('LineMaps/linesum-N_II-6583.fits')[0]
    x, y, w = hdu_true.data, hdu_syn.data, hduf.data
    xmin, xmax = ymin, ymax = pars['min'], pars['max']
    # mask out silly values
    m = np.isfinite(x) & np.isfinite(y/x) & (np.abs(np.log10(y/x)) < 1.0)
    if regionfile is not None:
        m = m & pyregion.open(regionfile).get_mask(hdu=hduf)
    
    H, xedges, yedges = np.histogram2d(x[m], y[m], 50,
                                       [[xmin, xmax], [ymin, ymax]],
                                       weights=w[m])
    ratiotext = '{} / {}'.format(pars['line1'], pars['line2'])
    fig, ax = plt.subplots(1, 1)
    ax.imshow((H.T)**(1.0/GAMMA), extent=[xmin, xmax, ymin, ymax],
              interpolation='nearest', aspect='auto', origin='lower', 
              cmap=cmap, alpha=1.0)
    ax.plot([xmin, xmax], [xmin, xmax], '-', alpha=1.0,
            lw=1, c='r', label=None)
    ax.set_xlabel('MUSE spectrum-derived line ratio')
    ax.set_ylabel('MUSE synthetic WFC3 filter-derived line ratio')
    ax.set_xlim(xmin, xmax)
    ax.set_ylim(ymin, ymax)
    ax.text(0.5, 0.15, ratiotext,
            horizontalalignment='center', transform=ax.transAxes)
    ax.text(0.5, 0.05, titles_from_extra[extra],
            horizontalalignment='center', transform=ax.transAxes)

    fig.set_size_inches(4.5, 4.5)
    fig.tight_layout(pad=2)
    fig.savefig(pltname)

    return pltname


if __name__ == '__main__':
    try:
        regionfile = sys.argv[1]
    except:
        regionfile = None
    for ratio_name, pars in plotpars.items():
        for extra in '', '-naive', '-flat':
            print(histogram_ratio_images(ratio_name, pars, extra, regionfile))
