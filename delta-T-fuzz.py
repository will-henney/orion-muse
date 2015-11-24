from __future__ import print_function
import sys
import glob
import numpy as np
from astropy.io import fits

prefix = 'Linemaps/'

try:
    quantity = sys.argv[1]
except:
    sys.exit('Usage: {} (Te|Te-iii|Ne|Ne-iii)'.format(sys.argv[0]))

obs_pattern = prefix + 'muse-derived-{}-bin???.fits'.format(quantity)
obs_datafiles = glob.glob(obs_pattern)

for obs_fn in obs_datafiles:
    obs_hdu = fits.open(obs_fn)['SCALED'].data
    fuzz_pattern = obs_fn.replace('-bin', '-fuzz???-bin')
    fuzz_datafiles = glob.glob(fuzz_pattern)
    for fuzz_fn in fuzz_datafiles:
        fuzz_hdu = fits.open(fuzz_fn)['SCALED']
        if 'Te' in fuzz_fn:
            # calculate (T - T0)/T0
            fuzz_hdu.data = (fuzz_hdu.data - obs_hdu.data)
        else:
            # calculate log10(N/N0)
            fuzz_hdu.data = np.log10(fuzz_hdu.data/obs_hdu.data)

        delta_fn = fuzz_fn.replace('derived-Te', 'delta-Te').replace('derived-Ne', 'log10-n')
        print('Writing', delta_fn)
        fuzz_hdu.writeto(delta_fn, clobber=True)
