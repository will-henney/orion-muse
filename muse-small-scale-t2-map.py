from __future__ import print_function
import sys
import numpy as np
from rebin_utils import oversample
from astropy.io import fits

nlist = [4, 8, 16, 32, 64, 128]
for suffix in '', '-iii':
    fn = 'LineMaps/muse-derived-Te{}-bin001.fits'.format(suffix)
    hdu = fits.open(fn)['SCALED']
    ny, nx = hdu.data.shape
    print(ny, nx)
    for nbin in nlist:
        sbin = 'bin{:03d}'.format(nbin)
        mfn = fn.replace('bin001', sbin)
        T0 = fits.open(mfn)['SCALED'].data
        efn = fn.replace('bin', 'STD-bin')
        std = fits.open(efn)['SCALED'].data
        dt = (hdu.data - T0)/T0
        edt = std/T0
        dt4d = dt.reshape((nbin, ny//nbin, nbin, nx//nbin))
        t2 = np.nanmedian(dt4d**2, axis=(0, 2))
        t2 = oversample(t2, nbin)
        edt4d = edt.reshape((nbin, ny//nbin, nbin, nx//nbin))
        et2 = np.nanmedian(edt4d**2, axis=(0, 2))
        et2 = oversample(et2, nbin)
        t2 -= et2
        # sfn = mfn.replace('bin', 'SN-bin')
        # sn = fits.open(sfn)['SCALED'].data
        # t2 -= 1./sn**2
        out_fn = mfn.replace('derived-Te', 'finescale-t2')
        print('Writing', out_fn)
        fits.PrimaryHDU(header=hdu.header, data=t2).writeto(out_fn,
                                                            clobber=True)
        out_fn = mfn.replace('derived-Te', 'finescale-et2')
        print('Writing', out_fn)
        fits.PrimaryHDU(header=hdu.header, data=et2).writeto(out_fn,
                                                            clobber=True)
        out_fn = mfn.replace('derived-Te', 'finescale-dt')
        print('Writing', out_fn)
        fits.PrimaryHDU(header=hdu.header, data=dt).writeto(out_fn,
                                                            clobber=True)
