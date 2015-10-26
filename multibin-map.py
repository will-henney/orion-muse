from __future__ import print_function
import sys
import numpy as np
sys.path.append('/Users/will/Work/RubinWFC3/Tsquared')
from rebin_utils import downsample, oversample
from astropy.io import fits

nlist = [1, 2, 4, 8, 16, 32, 64, 128, 256]
mingoods = [2, 2, 2, 2, 2, 2, 2, 2, 2]

try: 
    infile = sys.argv[1]
except:
    sys.exit('Usage: {} FITSFILE'.format(sys.argv[0]))


hdu = fits.open(infile)[0]
hdr = hdu.header
im = hdu.data
w = np.ones_like(hdu.data)

continuum = fits.open('muse-hr-image-wfc3-f547m.fits')['DATA'].data
starmask = continuum > 30
m =  np.isfinite(hdu.data) & (hdu.data > 0.0) & (~starmask)
for n, mingood in zip(nlist, mingoods):
    im[~m] = 0.0
    outfile = infile.replace('.fits', '-bin{:03d}.fits'.format(n))
    print('Saving', outfile)
    # Save both the scaled image and the weights, but at the full resolution
    fits.HDUList([
        fits.PrimaryHDU(),
        fits.ImageHDU(data=oversample(im, n), header=hdr, name='scaled'),
        fits.ImageHDU(data=oversample(w, n), header=hdr, name='weight'),
    ]).writeto(outfile, clobber=True)
    # Now do the rebinning by a factor of two
    [im,], m, w = downsample([im,], m, weights=w, mingood=mingood)
