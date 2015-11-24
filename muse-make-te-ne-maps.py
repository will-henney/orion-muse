import sys
import numpy as np
from astropy.io import fits
from derive_ne_te_1phase import T_den_from_rsii_rnii

prefix = 'Linemaps/'
hduA = fits.open(prefix + "ratio-6716-6731.fits")[0]
hduB = fits.open(prefix + "ratio-5755-6583-deredden-2874.fits")[0]

Te = np.empty_like(hduA.data)
Ne = np.empty_like(hduA.data)
m = np.isfinite(hduA.data) & np.isfinite(hduB.data) & (hduA.data > 0) & (hduB.data > 0)
Te[m], Ne[m] = T_den_from_rsii_rnii(hduA.data[m], hduB.data[m])
Te[~m], Ne[~m] = np.nan, np.nan
fits.PrimaryHDU(header=hduA.header, data=Te).writeto(
    prefix + 'muse-derived-Te.fits', clobber=True)
fits.PrimaryHDU(header=hduA.header, data=Ne).writeto(
    prefix + 'muse-derived-Ne.fits', clobber=True)
