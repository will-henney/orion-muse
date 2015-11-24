from __future__ import print_function
import sys
import glob
from distutils.dep_util import newer, newer_group
import numpy as np
from skimage.morphology import binary_dilation, binary_erosion, square, diamond
from skimage.filters.rank import modal
from astropy.io import fits
from rebin_utils import oversample

def minify(a, n):
    return a[::n, ::n]

try: 
    fileroot = sys.argv[1]
    target_signal_to_noise = int(sys.argv[2])
except IndexError:
    sys.exit('Usage: {} FILEROOT S/N')

prefix = 'LineMaps/'
nlist = [1, 2, 4, 8, 16, 32, 64, 128, 256]
# selection element for filtering out isolated pixels
element = square(3)
for n in reversed(nlist):
    fn = prefix + '{}-SN-bin{:03d}.fits'.format(fileroot, n)
    out_fn = fn.replace('SN', 'mask-SN{:04d}'.format(target_signal_to_noise))
    if not newer(fn, out_fn):
        print(out_fn, 'already up to date - skipping')
        continue
    print('Extracting mask from', fn)
    try:
        hdu = fits.open(fn)['SCALED']
    except IOError:
        sys.exit(fn + ' not found. Try running multibin-signal-to-noise.py first.')
    hdr = hdu.header
    signal_to_noise = hdu.data
    mask = signal_to_noise >= target_signal_to_noise
    # Shrink down to its true size at the binning we are at
    mask = minify(mask, n).astype(np.uint8)
    # Eliminate small islands
    mask = mask & modal(mask, element)
    # mask = binary_erosion(mask, element)
    # expand borders slightly
    # mask = binary_dilation(mask, element)

    # Expand back to the size of the full-res image
    mask = oversample(mask, n).astype(np.uint8)
    # Save to FITS file
    fits.PrimaryHDU(header=hdr, data=mask).writeto(out_fn, clobber=True)
