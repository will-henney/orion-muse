from __future__ import print_function
import sys
import numpy as np
import pyneb
from astropy.io import fits

# Set up Blagrave 2007 extinction law
REDCORR = pyneb.extinction.red_corr.RedCorr(
    law='CCM89 Bal07', R_V=5.5, cHbeta=1.0)

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
        lineid = sys.argv[1]
        ion_string, wav_string = lineid.split('-')
    except IndexError:
        print('Usage: {} LINEID'.format(sys.argv[0]))
    try:
        suffix = sys.argv[2] + '.fits'
        iHDU = 'SCALED'
    except:
        suffix = '.fits'
        iHDU = 0

    hb_ha = fits.open('Linemaps/ratio-4861-6563' + suffix)[iHDU].data
    ha_h9229 = fits.open('Linemaps/ratio-6563-9229' + suffix)[iHDU].data 
    chb = CHb_from_RHbHa(hb_ha)
    chb2 = CHb_from_R6563_9229(ha_h9229)

    wav = int(wav_string)
    clam = (1.0 + flambda(wav))*chb
    clam2 = (1.0 + flambda(wav))*chb2

    for prefix in 'continuum', 'linesum':
        fn = 'Linemaps/{}-{}'.format(prefix, lineid) + suffix
        try:
            hdu = fits.open(fn)[iHDU]
        except IOError:
            print('Skipping', fn, iHDU)
            continue
        for c, newsuff in [[clam, '-excorr'],
                           [clam2, '-excorr2']]:
            fn_new = fn.replace(suffix, newsuff + suffix)
            fits.PrimaryHDU(
                data=hdu.data*10**c,
                header=hdu.header).writeto(fn_new, clobber=True)
            print(fn_new)
