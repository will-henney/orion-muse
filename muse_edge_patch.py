import sys
from astropy.io import fits
from astropy import wcs
from astropy import units as u
import numpy as np

sections = np.arange(8, dtype=int)
NV = 702
# Extract cube that is +/- dwav around the wavsec boundary 
dwav = 12.0

try:
    quantity = sys.argv[1]
    section = int(sys.argv[2])
    assert quantity in ['data', 'variance']
    # makes no sense for 0th section
    assert section in sections and section > 0
except (IndexError, AssertionError) as e:
    sys.exit('Usage: {} (data|variance) SECTION'.format(sys.argv[0]))


if quantity == 'variance':
    cubename = 'STAT'
else:
    cubename = 'DATA'

hdulist = fits.open('DATA/DATACUBEFINALuser_20140216T010259_78380e1d.fits')
cube = hdulist[cubename]
w = wcs.WCS(cube.header)
wav0 = cube.header['CRVAL3'] + cube.header['CD3_3']*NV*section
wavs = (wav0 + dwav*np.array([-1.0, 1.0]))*u.Angstrom.to(u.m)
print('Wavelength window', wavs)
_, _, pixels = w.wcs_world2pix([0, 0], [0, 0], wavs, 0)
k1, k2 = int(pixels[0]), int(pixels[1]) + 2
window = slice(k1, k2), slice(None), slice(None)
fn = 'muse-hr-{}-wavsec-edge{}{}.fits'.format(quantity, section - 1, section)
hdr = cube.header.copy()
hdr.update(w.slice(window).to_header())

print('Writing', fn)
fits.PrimaryHDU(header=hdr, data=cube.data[window]).writeto(fn, clobber=True)
