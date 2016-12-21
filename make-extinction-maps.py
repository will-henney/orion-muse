from __future__ import print_function
import sys
import numpy as np
import pyneb
from astropy.io import fits

# Set up Blagrave 2007 extinction law
REDCORR = pyneb.extinction.red_corr.RedCorr(
    law='CCM89 Bal07', R_V=5.5, cHbeta=1.0)

WAV_PASCHEN_BETA = 12820 

def flambda(wav):
    """Find [(A_lam / A_Hb) - 1] as function of wavelength `wav`

    This is the same as given in Table 2 of Blagrave et al (2007)

    """
    return np.log10(REDCORR.getCorrHb(wav))


def CHb_from_RHbHa(RHbHa, balmer0=2.874):
    """Find base-10 extinction at H beta from balmer decrement `RHbHa`

    Assumes that the intrinsic Balmer decrement is `balmer0`

    """
    return np.log10(balmer0*RHbHa) / flambda(6563)


def CHb_from_R6563_9229(RBaPa, RBaPa0=112.0):
    """Find base-10 extinction at H beta from 6563/9229 decrement `RBaPa`

    Assumes that the intrinsic decrement is `RBaPa0`

    """
    return np.log10(RBaPa/RBaPa0) / (flambda(9229) - flambda(6563))


if __name__ == '__main__':

    try:
        suffix = sys.argv[1] + '.fits'
        iHDU = 'SCALED'
    except:
        suffix = '.fits'
        iHDU = 0

    hdu = fits.open('Linemaps/ratio-4861-6563' + suffix)[iHDU]
    hb_ha = hdu.data
    ha_h9229 = fits.open('Linemaps/ratio-6563-9229' + suffix)[iHDU].data 
    chb = CHb_from_RHbHa(hb_ha)
    chb2 = CHb_from_R6563_9229(ha_h9229)

    # Find C(Pa b) / C(H b)
    paschen_b_factor = (1.0 + flambda(WAV_PASCHEN_BETA))
    print('C(Pa b) / C(H b) =', paschen_b_factor)

    for name, data in [['C_H_beta', chb],
                       ['C_Pa_beta', paschen_b_factor*chb],
                       ['C_H_beta_9229', chb2],
                       ['C_Pa_beta_9229', paschen_b_factor*chb2]]:
        fn = 'Linemaps/' + name + suffix
        fits.PrimaryHDU(data=data, header=hdu.header).writeto(fn, clobber=True)
        print(fn)
