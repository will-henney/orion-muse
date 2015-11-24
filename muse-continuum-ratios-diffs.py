import numpy as np
from astropy.io import fits

def cfile(wav, binning):
    return 'LineMaps/continuum-average-{}-bin{:03d}.fits'.format(wav, binning)

def rfile(wav1, wav2, binning):
    return 'LineMaps/continuum-ratio-{}-{}-bin{:03d}.fits'.format(wav1, wav2, binning)

def dfile(wav1, wav2, binning):
    return 'LineMaps/continuum-diff-{}-{}-bin{:03d}.fits'.format(wav1, wav2, binning)


pairs = [
    ['4000', '6000'],
    ['4000', '5000'],
    ['5000', '6000'],
    ['6000', '7000'],
    ['7000', '8000'],
    ['8000', '9000'],
    ]
binnings = [1, 2, 4, 8, 16, 32, 64, 128, 256]
for binning in binnings:
    for wavA, wavB in pairs:
        hduA = fits.open(cfile(wavA, binning))['SCALED']
        hduB = fits.open(cfile(wavB, binning))['SCALED']

        fn = rfile(wavA, wavB, binning)
        print('Writing', fn)
        fits.PrimaryHDU(header=hduA.header,
                        data=hduA.data/hduB.data).writeto(fn, clobber=True)

        fn = dfile(wavA, wavB, binning)
        print('Writing', fn)
        fits.PrimaryHDU(header=hduA.header,
                        data=hduA.data-hduB.data).writeto(fn, clobber=True)
