from __future__ import print_function
import sys
import glob
import numpy as np
from astropy.io import fits
try: 
    fileroot = sys.argv[1]
    target_signal_to_noise = int(sys.argv[2])
    maskroot = sys.argv[3]
    newsuff = sys.argv[4]
except IndexError:
    sys.exit('Usage: {} FILEROOT S/N MASKROOT NEWSUFF')

prefix = 'LineMaps/'
snlabel = 'SN{:04d}'.format(target_signal_to_noise)
outim = None
nlist = [1, 2, 4, 8, 16, 32, 64, 128, 256]
for n in reversed(nlist):
    # Read in data
    fn = prefix + '{}-bin{:03d}.fits'.format(fileroot, n)
    hdu = fits.open(fn)['SCALED']
    if outim is None:
        # One-time setup to do on first iteration
        hdr = hdu.header
        # set up image for output
        outim = np.empty_like(hdu.data)

    # Read in mask
    mfn = fn.replace(fileroot,
                     '{}-mask-{}'.format(maskroot, snlabel))
    try:
        mask = fits.open(mfn)['SCALED'].data.astype(np.bool)
    except IOError:
        sys.exit(mfn + ' not found. Try running multibin-mask-s-n.py first.')
    # mask = mask & (hdu.data > 0.0)

    # Paste into image and into saved s/n
    outim[mask] = hdu.data[mask]

out_fn = '{}{}-multibin-{}.fits'.format(prefix, fileroot, newsuff)
fits.PrimaryHDU(header=hdr, data=outim).writeto(out_fn, clobber=True)
print(out_fn)
