from __future__ import print_function
from misc_utils import sanitize_string
from astropy.io import fits
import nebulio
print(nebulio.__version__)

# Mean and std of the color terms ktwiddle
# This is copied from ratio-sensitivity-nebulio.py in the Tsquared folder
COLOR_TERMS_MEAN_SIG = {
    "FQ674N": (1.00, 0.10), 
    "F673N":  (1.00, 0.10), 
    "FQ672N": (0.99, 0.10), 
    "F658N":  (1.04, 0.11), 
    "F656N":  (1.04, 0.11), 
    "FQ575N": (0.95, 0.03), 
    "F547M":  (1.03, 0.01), 
    "F502N":  (1.10, -0.04), 
    "F487N":  (1.06, -0.06), 
#    "F487N":  (1.65, -0.06), 
    "FQ437N": (1.41, -0.13), 
    "FQ436N": (1.87, -0.08), 
}
def get_color_term(fname="FQ575N", kshift=0.0):
    """Return color term for filter `fname`, shifted `kshift` stdevs from mean"""
    mean, sigma = COLOR_TERMS_MEAN_SIG[fname]
    return mean + kshift*sigma

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
    hdu = fits.open('muse-hr-image-wfc3-{}.fits'.format(fn.lower()))['DATA']
    return hdu.data

def get_fits_header(fn='FQ575N'):
    hdu = fits.open('muse-hr-image-wfc3-{}.fits'.format(fn.lower()))['DATA']
    return hdu.header

heliocentric_correction = -16.2
default_velocity =  25.0 + heliocentric_correction
default_width = 20.0

for ratio_name, filterset in filtersets.items():
    FI, FII, FIII = [filterset[J] for J in ("I", "II", "III")]
    print(FI, FII, FIII)
    RI = get_fits_data(FI)
    RII = get_fits_data(FII)
    RIII = get_fits_data(FIII)
    print(RIII[::100, ::100])
    lineids = filterset['line1'], filterset['line2']
    bpnames = ['wfc3,uvis1,' + F for F in (FI, FII, FIII)]
    fset = nebulio.Filterset(bpnames, lineids,
                             velocity=default_velocity, fwhm_kms=default_width)
    kI = get_color_term(FI)/get_color_term(FIII)
    kII = get_color_term(FII)/get_color_term(FIII)
    print('Color terms:', kI ,kII)
    ratio = fset.find_line_ratio(rates=[RI, RII, RIII], colors=(kI, kII), naive=False)
    ratio_naive = fset.find_line_ratio(rates=[RI, RII, RIII], colors=(kI, kII), naive=True)
    # And if we just ignore the color terms
    ratio_flat = fset.find_line_ratio(rates=[RI, RII, RIII], colors=(1.0, 1.0), naive=False)

    outhdu = fits.PrimaryHDU(header=get_fits_header(), data=ratio)
    outhdu.writeto(
        'NebulioMUSE/synthetic-ratio-{}.fits'.format(ratio_name),
        clobber=True)

    outhdu = fits.PrimaryHDU(header=get_fits_header(), data=ratio_naive)
    outhdu.writeto(
        'NebulioMUSE/synthetic-naive-ratio-{}.fits'.format(ratio_name),
        clobber=True)

    outhdu = fits.PrimaryHDU(header=get_fits_header(), data=ratio_flat)
    outhdu.writeto(
        'NebulioMUSE/synthetic-flat-ratio-{}.fits'.format(ratio_name),
        clobber=True)

    # Take ratio of ratios between simulated and true
    # Ignore anything after the second dash in finding the true ratio name
    ratio_name_true = '-'.join(ratio_name.split('-')[:2])
    fn_true = 'LineMaps/ratio-{}.fits'.format(ratio_name_true)
    ratio_true = fits.open(fn_true)[0].data
    ratio_of_ratios = ratio/ratio_true
    outhdu = fits.PrimaryHDU(header=get_fits_header(), data=ratio_of_ratios)
    outhdu.writeto(
        'NebulioMUSE/synthetic-over-true-ratio-{}.fits'.format(ratio_name),
        clobber=True)
    # And the same for the naive/true
    ratio_of_ratios = ratio_naive/ratio_true
    outhdu = fits.PrimaryHDU(header=get_fits_header(), data=ratio_of_ratios)
    outhdu.writeto(
        'NebulioMUSE/synthetic-naive-over-true-ratio-{}.fits'.format(ratio_name),
        clobber=True)
    # And the same for the flat/true
    ratio_of_ratios = ratio_flat/ratio_true
    outhdu = fits.PrimaryHDU(header=get_fits_header(), data=ratio_of_ratios)
    outhdu.writeto(
        'NebulioMUSE/synthetic-flat-over-true-ratio-{}.fits'.format(ratio_name),
        clobber=True)
