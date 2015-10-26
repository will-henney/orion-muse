import sys
from astropy.io import fits

sys.path.append('/Users/will/Work/RubinWFC3/Tsquared')
from deredden import deredden_nii_ratio

hb_ha = fits.open('Linemaps/ratio-4861-6563.fits')[0].data
nii_hdu = fits.open('Linemaps/ratio-5755-6583.fits')[0]
nii_hdu.data = deredden_nii_ratio(nii_hdu.data, hb_ha)
nii_hdu.writeto('Linemaps/ratio-5755-6583-deredden-2874.fits', clobber=True)

hb_ha = fits.open('NebulioMUSE/synthetic-ratio-4861-6563.fits')[0].data
nii_hdu = fits.open('NebulioMUSE/synthetic-ratio-5755-6583.fits')[0]
nii_hdu.data = deredden_nii_ratio(nii_hdu.data, hb_ha)
nii_hdu.writeto('NebulioMUSE/synthetic-ratio-5755-6583-deredden-2874.fits', clobber=True)
