from __future__ import print_function
import sys
import glob
from distutils.dep_util import newer, newer_group
import numpy as np
from astropy.io import fits
try: 
    fileroot = sys.argv[1]
except IndexError:
    sys.exit('Usage: {} FILEROOT')

prefix = 'LineMaps/'
nlist = [1, 2, 4, 8, 16, 32, 64, 128, 256]
for n in reversed(nlist):
    fn = prefix + '{}-bin{:03d}.fits'.format(fileroot, n)
    out_fn = fn.replace('-bin', '-SN-bin')
    if not newer(fn, out_fn):
        print(out_fn, 'already up to date - skipping')
        continue
    hdu = fits.open(fn)['SCALED']
    hdr = hdu.header
    fuzz_pattern = fn.replace('-bin', '-fuzz???-bin')
    fuzz_files = glob.glob(fuzz_pattern)
    if len(fuzz_files) == 0:
        sys.exit('No fuzzed files found. '
                 'Try running extract-em-line-fuzz.py and all-fuzz-multibin.sh ')
    hdr['NFUZZ'] = len(fuzz_files), 'Number of fuzzed images used to estimate noise'
    fuzz_deltas = []
    for fuzz_fn in fuzz_files:
        fuzz_hdu = fits.open(fuzz_fn)['SCALED']
        fuzz_deltas.append(fuzz_hdu.data - hdu.data)
    fuzz_std = np.sqrt(np.mean(np.square(np.dstack(fuzz_deltas)), axis=-1))
    signal_to_noise = np.abs(hdu.data)/fuzz_std
    # Write out signal/noise
    fits.PrimaryHDU(header=hdr, data=signal_to_noise).writeto(out_fn, clobber=True)
    # Write out noise too
    out_fn = fn.replace('-bin', '-STD-bin')
    fits.PrimaryHDU(header=hdr, data=fuzz_std).writeto(out_fn, clobber=True)
