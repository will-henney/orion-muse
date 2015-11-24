from __future__ import print_function
import sys
import numpy as np
from astropy.io import fits
from matplotlib import pyplot as plt
from wcsaxes import WCS
# import seaborn as sns

def fig_ax_im_from_fits(fn, vmin=None, vmax=None, cmap=plt.cm.gray):
    '''Make a basic grayscale plot of a FITS file'''
    #sns.set_style("white")

    # Try first HDU
    hdu = fits.open(fn)[0]
    if hdu.data is None:
        # Or failing that, 2nd HDU
        hdu = fits.open(fn)[1]
    wcs = WCS(hdu.header)
    fig = plt.figure()
    ax = fig.add_axes([0.18, 0.1, 0.8, 0.8], projection=wcs)
    tr = ax.get_transform(wcs)
    ra, dec = ax.coords
    dec.set_major_formatter('dd:mm:ss.s')
    ra.set_major_formatter('hh:mm:ss.ss')
    if vmax is None:
        vmax = hdu.data.max()
    if vmin is None:
        vmin = 0.0
    im = ax.imshow(hdu.data, vmin=vmin, vmax=vmax, cmap=cmap,
                   origin='lower', interpolation='nearest')
    plt.colorbar(im)
    ax.coords.grid(color='red', alpha=0.2, lw=0.2, linestyle='solid')
    ax.set_xlabel('RA, J2000')
    ax.set_ylabel('Dec, J2000')
    ax.set_title(fn, fontsize='small')

    return fig, ax, im


if __name__ == '__main__':
    try:
        fn = 'LineMaps/' + sys.argv[1] + '.fits'
    except IndexError:
        fn = 'LineMaps/linesum-Fe_III-4658-multibin-SN0010.fits'

    try:
        vmin, vmax = float(sys.argv[2]), float(sys.argv[3])
    except IndexError:
        vmin, vmax = None, None
        
    fig, ax, im = fig_ax_im_from_fits(fn, vmin, vmax)
    figfile = fn.replace('.fits', '.jpg')
    fig.set_size_inches(8, 6)
    fig.savefig(figfile, dpi=200)
    print(figfile)
