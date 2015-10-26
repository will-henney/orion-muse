import sys
import numpy as np
from astropy.io import fits
t2dir = '/Users/will/Work/RubinWFC3/Tsquared'
sys.path.append(t2dir)
from derive_ne_te_1phase import T_den_from_rcliii_rsiii

prefix = 'Linemaps/'
nlist = [1, 2, 4, 8, 16, 32, 64, 128, 256]

for n in nlist:
    suffix = '-bin{:03d}'.format(n)

    hduA = fits.open(prefix + "ratio-5538-5518{}.fits".format(suffix))[1]
    hduB = fits.open(prefix + "ratio-6312-9069-deredden{}.fits".format(suffix))[1]

    Te = np.empty_like(hduA.data)
    Ne = np.empty_like(hduA.data)
    m = (np.isfinite(hduA.data) & np.isfinite(hduB.data)
         & (hduA.data > 0) & (hduB.data > 0))
    Te[m], Ne[m] = T_den_from_rcliii_rsiii(hduA.data[m], hduB.data[m])
    Te[~m], Ne[~m] = np.nan, np.nan
    fits.PrimaryHDU(header=hduA.header, data=Te).writeto(
        prefix + 'muse-derived-Te-iii{}.fits'.format(suffix), clobber=True)
    fits.PrimaryHDU(header=hduA.header, data=Ne).writeto(
        prefix + 'muse-derived-Ne-iii{}.fits'.format(suffix), clobber=True)
