from __future__ import print_function
import sys
import glob
import numpy as np
from skimage.morphology import erosion, dilation
from astropy.io import fits
try: 
    fileroot = sys.argv[1]
    target_signal_to_noise = int(sys.argv[2])
except IndexError:
    sys.exit('Usage: {} FILEROOT S/N')

prefix = 'LineMaps/'
snlabel = 'SN{:04d}'.format(target_signal_to_noise)
outim = None
nlist = [1, 2, 4, 8, 16, 32, 64, 128, 256]
for n in reversed(nlist):
    # Read in data
    fn = prefix + '{}-bin{:03d}.fits'.format(fileroot, n)
    hdu = fits.open(fn)['SCALED']
    # Read in s/n too
    sn_hdu = fits.open(fn.replace('-bin', '-SN-bin'))['SCALED']
    if outim is None:
        # One-time setup to do on first iteration
        hdr = hdu.header
        # set up image for output
        outim = np.empty_like(hdu.data)
        # and for s/n achieved
        out_sn = np.empty_like(hdu.data)

    # Read in mask
    mfn = fn.replace('-bin', '-mask-{}-bin'.format(snlabel))
    try:
        mask = fits.open(mfn)['SCALED'].data.astype(np.bool)
    except IOError:
        sys.exit(mfn + ' not found. Try running multibin-mask-s-n.py first.')
    # mask = mask & (hdu.data > 0.0)

    # Paste into image and into saved s/n
    outim[mask] = hdu.data[mask]
    out_sn[mask] = sn_hdu.data[mask]

out_fn = '{}{}-multibin-{}.fits'.format(prefix, fileroot, snlabel)
fits.PrimaryHDU(header=hdr, data=outim).writeto(out_fn, clobber=True)
fits.PrimaryHDU(header=hdr, data=out_sn).writeto(
    out_fn.replace('.fits', '-SN.fits'), clobber=True)
