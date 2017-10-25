from __future__ import print_function
import sys
import numpy as np
from astropy.table import Table
from astropy.io import fits
from astropy.wcs import WCS
from misc_utils import sanitize_string
from extract_utils import (find_wavsec, trim_to_window,
                           extract_line_maps, get_masks, linetab)
from fuzz_utils import fuzz


def save_linemap_files(
        wav, species,
        mapdir='LineMaps', usecont=[1, 1], specmask=None):
    full_width = 24.0  # Angstrom
    wavsec = find_wavsec(wav)
    wavsec_blue = find_wavsec(wav - full_width/2)
    wavsec_red = find_wavsec(wav + full_width/2)
    if (wavsec_blue != wavsec) or (wavsec_red != wavsec):
        print('Uh, oh - line straddles wavsecs', wavsec_blue, wavsec, wavsec_red)
        wavsecs = set([wavsec_blue, wavsec, wavsec_red])
        fn = 'muse-hr-data-wavsec-edge{}{}.fits'.format(*wavsecs)
    else:
        fn = 'muse-hr-data-wavsec{}.fits'.format(wavsec)

    try:
        hdulist = fits.open(fn)
    except:
        # Maybe we have the cube in the BigFiles folder
        hdulist = fits.open('BigFiles/' + fn)
    hdu = hdulist['DATA']
    wfull = WCS(hdu.header)

    # Read in variance for fuzzing
    vfn = fn.replace('data', 'variance')
    try: 
        variance = fits.open(vfn)['STAT'].data
    except:
        variance = fits.open('BigFiles/' + vfn)['STAT'].data

    cube, w = trim_to_window(wav, hdu.data, wfull, dwav=full_width/2)
    varcube, _ = trim_to_window(wav, variance, wfull, dwav=full_width/2)

    fuzzies = fuzz(cube, varcube, nfuzz)

    for ifuzz, fuzzcube in enumerate(fuzzies):
        maps, spec =  extract_line_maps(wav, fuzzcube, w, usecont, specmask=specmask)
        wavid = str(int(wav+0.5))
        fuzzid = '{:03d}'.format(ifuzz)
        # Save the maps to FITS file 
        for mapid, mapdata in maps.items():
            mhdu = fits.PrimaryHDU(header=w.celestial.to_header(), data=mapdata)
            mname = '{}/{}-{}-{}-fuzz{}.fits'.format(mapdir, mapid,
                                                     species, wavid, fuzzid)
            mhdu.writeto(mname, clobber=True)
        sname = '{}/spec1d-{}-{}-fuzz{}.tab'.format(mapdir, species, wavid, fuzzid)
        # And save the spectrum to TSV file
        Table(spec).write(sname, format='ascii.tab')


try:
    # Optionally specify a single line
    wav_wanted = int(sys.argv[1])
except:
    wav_wanted = None

try:
    # Optionally specify number of fuzzies
    nfuzz = int(sys.argv[2])
except:
    nfuzz = 10

masks = get_masks()
for row in linetab:
    if wav_wanted is not None and wav_wanted != int(0.5 + row['wav0']):
        # jump all unwanted lines if command line argument was given
        continue                
    print(row['Ion'], row['wav0'])
    save_linemap_files(
        row['wav0'], sanitize_string(row['Ion']),
        usecont=[row['blue cont'], row['red cont']],
        specmask=masks[row['Class']],
    )
