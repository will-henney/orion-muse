import sys
import numpy as np
from astropy.io import fits
from derive_ne_te_1phase import T_den_from_rsii_rnii
from deredden import deredden_nii_ratio

prefix = 'LineMaps/'
nlist = [1, 2, 4, 8, 16, 32, 64, 128, 256]

try:
    extra = sys.argv[1]
except IndexError:
    extra = ''

for n in nlist:
    suffix = '{}-bin{:03d}'.format(extra, n)

    hduA = fits.open(prefix + "ratio-6716-6731{}.fits".format(suffix))[1]
    hduB = fits.open(prefix + "ratio-5755-6583{}.fits".format(suffix))[1]
    hb_ha = fits.open(prefix + "ratio-4861-6563{}.fits".format(suffix))[1].data
    hduB.data = deredden_nii_ratio(hduB.data, hb_ha)
    hduB.writeto(prefix + 'ratio-5755-6583{}-deredden-2874.fits'.format(suffix),
                 clobber=True)

    Te = np.empty_like(hduA.data)
    Ne = np.empty_like(hduA.data)
    m = (np.isfinite(hduA.data) & np.isfinite(hduB.data)
         & (hduA.data > 0) & (hduB.data > 0))
    Te[m], Ne[m] = T_den_from_rsii_rnii(hduA.data[m], hduB.data[m])
    Te[~m], Ne[~m] = np.nan, np.nan
    fits.PrimaryHDU(header=hduA.header, data=Te).writeto(
        prefix + 'muse-derived-Te{}.fits'.format(suffix), clobber=True)
    fits.PrimaryHDU(header=hduA.header, data=Ne).writeto(
        prefix + 'muse-derived-Ne{}.fits'.format(suffix), clobber=True)
