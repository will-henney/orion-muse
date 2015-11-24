import sys
import numpy as np
from astropy.io import fits

# Blagrave extinction law interpolated on table
# -0.138 + (6312 - 5876)*(-0.218 - (-0.138))/(6548 - 5876)
f_6312 = -0.1899
# -0.309 + (9069 - 7330)*(-0.507 - (-0.309))/(9229 - 7330)
f_9069 = -0.4903
f_6563 = -0.220
f_9229 = -0.507

F_siii_ha_pa9 = (f_6312 - f_9069) / (f_6563 - f_9229) # Should be 1.0467

def deredden_siii_ratio(rsiii, rha_pa9, ha_pa9_intrinsic=113.0):
    """Uses the Blagrave reddening law"""
    return rsiii*(ha_pa9_intrinsic/rha_pa9)**F_siii_ha_pa9

if __name__ == '__main__':
    try:
        suffix = '-' + sys.argv[1]
    except IndexError:
        suffix = ''
      
    rha_pa9 = fits.open('LineMaps/ratio-6563-9229{}.fits'.format(suffix))['SCALED'].data
    siii_hdu = fits.open('LineMaps/ratio-6312-9069{}.fits'.format(suffix))['SCALED']
    siii_hdu.data = deredden_siii_ratio(siii_hdu.data, rha_pa9)
    siii_hdu.writeto('LineMaps/ratio-6312-9069-deredden{}.fits'.format(suffix),
                    clobber=True)
