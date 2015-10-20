from __future__ import print_function
import sys
from astropy.io import fits
def save_line_ratio_map(line1, line2, mapdir='LineMaps'):
    fn1 = '{}/linesum-{}.fits'.format(mapdir, line1)
    fn2 = '{}/linesum-{}.fits'.format(mapdir, line2)
    hdu1 = fits.open(fn1)[0]
    hdu2 = fits.open(fn2)[0]
    wav1 = line1.split('-')[-1]
    wav2 = line2.split('-')[-1]
    hdu1.data /= hdu2.data
    hdu1.writeto('{}/ratio-{}-{}.fits'.format(mapdir, wav1, wav2))

if __name__ == '__main__':
    try:
        line1 = sys.argv[1]
        line2 = sys.argv[2]
    except IndexError:
        sys.exit('Usage: {} LINE1 LINE2'.format(sys.argv[0]))

    save_line_ratio_map(line1, line2)
