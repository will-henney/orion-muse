import numpy as np
from astropy.modeling import models, fitting
from astropy.io import fits
from astropy.wcs import WCS
import pyregion

ZMIN, ZMAX = 0.5, 1.5
def fit_surf_to_image(hdu, goodregion=None, badregion=None):
    wcs = WCS(hdu.header)
    ny, nx = hdu.data.shape

    xpix, ypix = np.meshgrid(np.arange(nx), np.arange(ny))
    x, y = wcs.all_pix2world(xpix, ypix, 0)
    z = hdu.data
    z[z < ZMIN] = np.nan
    z[z > ZMAX] = np.nan
    mask = np.isfinite(z)
    if goodregion is not None:
        mask = mask & goodregion.get_mask(hdu=hdu)
    if badregion is not None:
        mask = mask & (~badregion.get_mask(hdu=hdu))
      
    p_init = models.Polynomial2D(degree=1)
    fit_p = fitting.LinearLSQFitter()
    p = fit_p(p_init, x[mask], y[mask], z[mask])
    zfit = p(x, y)
    zfit_mean = np.nanmean(zfit)
    zfit /= zfit_mean
    print(p.parameters/zfit_mean)
    return zfit


    
FILTERS = ['FQ575N', 'FQ672N', 'FQ674N', 'F656N', 'F658N', 'F487N', 'F547M', 'F673N']

T2DIR = '/Users/will/Work/RubinWFC3/Tsquared/'

if __name__ == '__main__':
    for filter_ in FILTERS:
        fn = 'wfc3-over-muse-new-calib-ratio-{}.fits'.format(filter_)
        hdu, = fits.open(fn)
        try:
            goodregion = pyregion.open('{}-good.reg'.format(filter_))
        except:
            goodregion = None
        try:
            badregion = pyregion.open('{}-bad.reg'.format(filter_))
        except:
            badregion = None
        print(filter_)
        surface = fit_surf_to_image(hdu, goodregion, badregion)
        hdu.data /= surface
        hdu.writeto(fn.replace('-new-', '-correct-'), clobber=True)
        hdu.data = surface
        hdu.writeto(fn.replace('-new-', '-surface-'), clobber=True)
