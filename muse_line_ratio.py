from __future__ import print_function
import sys
from astropy.io import fits
def save_line_ratio_map(line1, line2, prefix='linesum', suffix='', mapdir='LineMaps'):
    fn1 = '{}/{}-{}{}.fits'.format(mapdir, prefix, line1, suffix)
    fn2 = '{}/{}-{}{}.fits'.format(mapdir, prefix, line2, suffix)
    print(fn1, fn2)
    hdu1 = fits.open(fn1)[0]
    if hdu1.data is None:
        hdu1 = fits.open(fn1)[1]
    hdu2 = fits.open(fn2)[0]
    if hdu2.data is None:
        hdu2 = fits.open(fn2)[1]
    wav1 = line1.split('-')[-1]
    wav2 = line2.split('-')[-1]
    hdu1.data /= hdu2.data
    if prefix == 'linesum':
        hdu1.writeto('{}/ratio-{}-{}{}.fits'.format(mapdir, wav1, wav2, suffix), clobber=True)
    else:
        hdu1.writeto('{}/ratio-{}-{}-{}.fits'.format(mapdir, prefix, wav1, wav2, suffix), clobber=True)

if __name__ == '__main__':
    try:
        line1 = sys.argv[1]
        line2 = sys.argv[2]
    except IndexError:
        sys.exit('Usage: {} LINE1 LINE2 [PREFIX] [SUFFIX]'.format(sys.argv[0]))

    try:
        prefix = sys.argv[3]
    except IndexError:
        prefix = 'linesum'

    try:
        suffix = '-' + sys.argv[4]
    except IndexError:
        suffix = ''

    save_line_ratio_map(line1, line2, prefix, suffix)
