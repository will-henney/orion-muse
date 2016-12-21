from __future__ import print_function
import sys
from misc_utils import sanitize_string
from astropy.io import fits
T2DIR = '/Users/will/Work/RubinWFC3/Tsquared/'
sys.path.append(T2DIR)
import nebulio
import nebulio_adjustments_muse
from nebulio_adjustments_muse import set_transmission
from new_color_terms import get_color_term

print(nebulio.__version__)

filtersets = {
    "5755-6583": {"line1": "[N II] 5755", "line2": "[N II] 6583",
                         "I": "FQ575N", "II": "F658N", "III": "F547M"},
    "6716-6731": {"line1": "[S II] 6716", "line2": "[S II] 6731",
                         "I": "FQ672N", "II": "FQ674N", "III": "F547M"}, 
    "6716-6731-N": {"line1": "[S II] 6716", "line2": "[S II] 6731",
                         "I": "FQ672N", "II": "FQ674N", "III": "F673N"}, 
    "4861-6563": {"line1": "H I 4861", "line2": "H I 6563",
                         "I": "F487N", "II": "F656N", "III": "F547M"},
}


def get_fits_data(fn='FQ575N'):
    hdu = fits.open('wfc3-new-resample-muse-{}.fits'.format(fn))[0]
    return hdu.data

def get_fits_header(fn='FQ575N'):
    hdu = fits.open('wfc3-new-resample-muse-{}.fits'.format(fn))[0]
    return hdu.header

heliocentric_correction = 22.40
default_velocity =  25.0 + heliocentric_correction
default_width = 20.0
all_of_the_filters = ['fq575n', 'fq672n', 'fq674n', 'f673n',
                      'f487n' , 'f656n' , 'f658n' , 'f547m']
CORRECT_LINEAR_GRADIENT = {
    'nominal': [],
    '2016-03': ['fq672n'],
    '2016-03-allgrads': all_of_the_filters,
    '2016-03-bias': all_of_the_filters,
}

for transmission_set in 'nominal', '2016-03-bias', '2016-03-allgrads':
    set_transmission(transmission_set=transmission_set)
    for ratio_name, filterset in filtersets.items():
        FI, FII, FIII = [filterset[J] for J in ("I", "II", "III")]
        RI = get_fits_data(FI)
        RII = get_fits_data(FII)
        RIII = get_fits_data(FIII)
        for F, R in zip([FI, FII, FIII], [RI, RII, RIII]):
            name3 = 'wfc3-over-muse-surface-calib-ratio-{}.fits'.format(F)
            if F.lower() in CORRECT_LINEAR_GRADIENT[transmission_set]:
                # Linear correction surface to fix flat field
                print('Gradient correction for', F, transmission_set)
                hdus = fits.open(name3)[0]
                R /= hdus.data

        lineids = filterset['line1'], filterset['line2']
        bpnames = ['wfc3,uvis1,' + F for F in (FI, FII, FIII)]
        fset = nebulio.Filterset(bpnames, lineids,
                                 velocity=default_velocity, fwhm_kms=default_width)
        kI = get_color_term(FI)/get_color_term(FIII)
        kII = get_color_term(FII)/get_color_term(FIII)
        ratio = fset.find_line_ratio(rates=[RI, RII, RIII],
                                     colors=(kI, kII), naive=False)
        outhdu = fits.PrimaryHDU(header=get_fits_header(), data=ratio)
        fn = 'NebulioMUSE/wfc3-resample-{}-{}.fits'.format(transmission_set, ratio_name)
        print(fn)
        outhdu.writeto(fn, clobber=True)
