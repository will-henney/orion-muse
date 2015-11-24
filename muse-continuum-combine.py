import glob
import numpy as np
from astropy.io import fits

wavrange = ['4?', '5?', '6?', '7?', '8[2-9]', '9?']
binnings = [1, 2, 4, 8, 16, 32, 64, 128, 256]

for binning in binnings:
    for wavpat in wavrange:
        wav000 = wavpat[0]
        pattern = 'LineMaps/continuum-*-{}??-bin{:03d}.fits'.format(wavpat, binning)
        cfiles = glob.glob(pattern)
        average = None
        for fn in cfiles:
            hdu = fits.open(fn)['SCALED']
            if average is None:
                average = np.zeros_like(hdu.data)
            average += hdu.data
        average /= len(cfiles)
        out_fn = 'LineMaps/continuum-average-{}000-bin{:03d}.fits'.format(wav000, binning)
        print('Saving', out_fn)
        fits.PrimaryHDU(header=hdu.header, data=average).writeto(out_fn, clobber=True)
