import numpy as np
from astropy.io import fits
def deredden_nii_ratio(rnii, rhbha, balmer0=2.874):
    """Uses the Blagrave reddening law"""
    chb = -np.log10(balmer0*rhbha) / 0.220
    return rnii*10**(0.099*chb)

if __name__ == '__main__':
    hb_ha = fits.open('newratio-4861-6563-F547M.fits')[0].data
    nii_hdu = fits.open('newratio-5755-6583-F547M.fits')[0]
    nii_hdu.data = deredden_nii_ratio(nii_hdu.data, hb_ha)
    nii_hdu.writeto('newratio-5755-6583-F547M-deredden-2874.fits', clobber=True)
