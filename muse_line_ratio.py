from __future__ import print_function
import sys
from astropy.io import fits
from distutils.dep_util import newer, newer_group

def save_line_ratio_map(line1, line2,
                        prefix='linesum', suffix='',
                        mapdir='LineMaps', force=False):
    fn1 = '{}/{}-{}{}.fits'.format(mapdir, prefix, line1, suffix)
    fn2 = '{}/{}-{}{}.fits'.format(mapdir, prefix, line2, suffix)
    wav1 = line1.split('-')[-1]
    wav2 = line2.split('-')[-1]
    if prefix == 'linesum':
        fn_out = '{}/ratio-{}-{}{}.fits'.format(
            mapdir, wav1, wav2, suffix)
    else:
        fn_out = '{}/ratio-{}-{}-{}{}.fits'.format(
            mapdir, prefix, wav1, wav2, suffix)

    if newer_group([fn1, fn2], fn_out) or force:
        print(fn1, fn2)
        hdu1 = fits.open(fn1)[0]
        if hdu1.data is None:
            hdu1 = fits.open(fn1)[1]
        hdu2 = fits.open(fn2)[0]
        if hdu2.data is None:
            hdu2 = fits.open(fn2)[1]
        hdu1.data /= hdu2.data
        hdu1.writeto(fn_out, clobber=True)
    else:
        print(fn_out, 'already up to date.  Use force=True to recreate regardless.')

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
