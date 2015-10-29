import sys
from astropy.io import fits
from astropy import wcs
import numpy as np

try:
    quantity = sys.argv[1]
    assert quantity in ['data', 'variance']
except (IndexError, AssertionError) as e:
    sys.exit('Usage: {} (data|variance)'.format(sys.argv[0]))


if quantity == 'variance':
    cubename = 'STAT'
else:
    cubename = 'DATA'
  
hdulist = fits.open('DATA/DATACUBEFINALuser_20140216T010259_78380e1d.fits')
cube = hdulist[cubename]
sections = np.arange(8, dtype=int)
NV = 702
k1_list = sections*NV
k2_list = k1_list + NV
wav0_list = cube.header['CRVAL3'] + cube.header['CD3_3']*NV*sections
for section, k1, k2, wav0 in zip(sections, k1_list, k2_list, wav0_list):
    fn = 'muse-hr-{}-wavsec{}.fits'.format(quantity, section)
    hdr = cube.header.copy()
    hdr['CRVAL3'] = wav0
    hdr['NAXIS3'] = NV
    print('Writing', fn)
    fits.PrimaryHDU(header=hdr, data=cube.data[k1:k2]).writeto(fn, clobber=True)
