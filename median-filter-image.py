import sys
import numpy as np
import skimage
import skimage.filters
import skimage.morphology
import skimage.exposure
from astropy.io import fits

try:
    fn = sys.argv[1]
except IndexError:
    sys.exit('Usage: {} FITSFILE [RADIUS]'.format(sys.argv[0]))

try:
    radius = float(sys.argv[2])
except IndexError:
    radius = 3.5

hdulist = fits.open(fn)
hdu = hdulist[0]
if hdu.data is None:
    # Fallback if first HDU has no image
    hdu = hdulist['SCALED']

mask = np.isfinite(hdu.data)
save_range = hdu.data[mask].min(), hdu.data[mask].max()
selem = skimage.morphology.disk(radius)
image = skimage.exposure.rescale_intensity(
    hdu.data,
    in_range=save_range,
    out_range=(0.0, 1.0))
image = skimage.filters.median(
    skimage.img_as_uint(image),
    selem=selem, mask=mask)
hdu.data = skimage.exposure.rescale_intensity(
    skimage.img_as_float(image),
    out_range=save_range)
suffix = '-median{:03d}.fits'.format(int(radius))
hdu.writeto(fn.replace('.fits', suffix), clobber=True)
