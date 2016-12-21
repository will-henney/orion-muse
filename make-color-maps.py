"""Generate E/W maps for all the non-target lines
"""
import numpy as np
from astropy.io import fits
import nebulio
from misc_utils import sanitize_string

def fix_up_lineid(lineid):
    """Convert '[O III] 5007' to 'O_III-5007', etc"""
    pieces = lineid.split(' ')
    wavid = pieces[-1]
    ionid = sanitize_string(' '.join(pieces[:-1]))
    return '-'.join([ionid, wavid])


target_line_id = "[N II] 5755"
non_target_line_ids = [
    'He I 5876',
    '[Fe III] 5270',
    '[Cl III] 5518',
    '[Cl III] 5538',
    '[N I] 5199',
    '[O III] 5007',
]
FILTER_ID = 'f547m'
SUFFIX = '-bin016'
velocity = 25.0 - 16.2
width = 20.0

sanid = fix_up_lineid(target_line_id)
fn_I ='LineMaps/continuum-{}{}.fits'.format(sanid, SUFFIX) 
contI_hdu = fits.open(fn_I)['SCALED']
fn_III = 'muse-hr-image-wfc3-{}{}.fits'.format(FILTER_ID, SUFFIX)
R_III_hdu = fits.open(fn_III)['SCALED']

bp = nebulio.Bandpass('wfc3,uvis1,' + FILTER_ID)

C_WFC3 = 0.0840241
em = nebulio.EmissionLine(target_line_id, velocity=velocity, fwhm_kms=width)
lamIlam_cont_I = contI_hdu.data * 1e-20*em.wav0[0] / (0.2/206265)**2
lamIlam_III = R_III_hdu.data / (C_WFC3*bp.Tm*bp.Wj)
ktwiddle = lamIlam_cont_I / lamIlam_III
kfile = 'Linemaps/ktwid-{}-{}{}.fits'.format(sanid, FILTER_ID, SUFFIX)
contI_hdu.data = ktwiddle
fits.PrimaryHDU(header=contI_hdu.header, data=ktwiddle).writeto(kfile, clobber=True)

one_plus_sum_E_over_W = np.ones_like(ktwiddle)
for line_id in [target_line_id] + non_target_line_ids:
    sanid = fix_up_lineid(line_id)
    ewfile = 'LineMaps/ew-{}{}.fits'.format(sanid, SUFFIX)
    em = nebulio.EmissionLine(line_id, velocity=velocity, fwhm_kms=width)
    W = bp.Wtwid(em)[0]
    Ehdu = fits.open(ewfile)['SCALED']
    Ehdu.data /= W
    Ehdu.writeto(ewfile.replace('ew-', 'E_over_W-'), clobber=True)
    if line_id != target_line_id:
        one_plus_sum_E_over_W += Ehdu.data

fits.PrimaryHDU(header=contI_hdu.header,
                data=one_plus_sum_E_over_W
).writeto(
    'LineMaps/one_plus_sum_E_over_W-{}{}.fits'.format(
        FILTER_ID, SUFFIX),
    clobber=True)

fits.PrimaryHDU(header=contI_hdu.header,
                data=ktwiddle*one_plus_sum_E_over_W
).writeto(kfile.replace('ktwid', 'knorm'), clobber=True)
