from __future__ import print_function
import sys
import numpy as np
from astropy.table import Table
from astropy.io import fits
from astropy.wcs import WCS
from astropy import units as u
from misc_utils import sanitize_string
from helio_utils import waves2vels


linetab = Table.read('basic-line-list.tab', format='ascii.tab')
wavsectab = Table.read('wavsec-startwavs.tab', format='ascii.tab')

def find_wavsec(wav):
    """Which wavsec chunk is this wavelength in"""
    return wavsectab['Section'][wav > wavsectab['CRVAL3']].max() 


def trim_to_window(wav, cube, wcube, dwav=12.0):
    """Trim a spectral `cube` with accompanying WCS `wcube` down to a
window of width 2*`dwav` Angstrom around a central wavelength `wav`.
Returns tuple of trimmed cube and trimmed WCS

    """
    wavs = (wav + dwav*np.array([-1, 1]))*u.Angstrom.to(u.m)
    _, _, [wavmin, wavmax] = wcube.wcs_pix2world([0, 0], [0, 0], [0, cube.shape[0]-1], 0)
    print('Requested wavelength window:', wavs)
    print('Available wavelength range:', wavmin, wavmax)
    _, _, pixels = wcube.wcs_world2pix([0, 0], [0, 0], wavs, 0)
    k1, k2 = int(pixels[0]), int(pixels[1])+2
    if k1 < 0:
        print('Adjusting blue edge from', k1, 'to', 0)
        k1 = 0
    window = slice(k1, k2), slice(None), slice(None)
    return cube[window], wcube.slice(window)


def extract_line_maps(wav, cube, wcube, helio_hdr,
                      usecont=[1, 1], dwav=4.0, dwavcont=8.0):
    """Extract line moments and continuum"""

    # outer range is wav +/- dwavcont
    wavs = (wav + dwavcont*np.array([-1, 1]))*u.Angstrom.to(u.m)
    _, _, (kout1, kout2) = wcube.all_world2pix([0, 0], [0, 0], wavs, 0)
    kout1, kout2 = int(kout1), int(kout2) + 2
    if kout1 < 0:
        print('Adjusting kout1 from', kout1, 'to 0')
        kout1 = 0

    # inner range is wav +/- dwav
    wavs = (wav + dwav*np.array([-1, 1]))*u.Angstrom.to(u.m)
    _, _, (kin1, kin2) = wcube.all_world2pix([0, 0], [0, 0], wavs, 0)
    kin1, kin2 = int(kin1), int(kin2) + 2
    print('kout1, kin1, kin2, kout2, nwav =', kout1, kin1, kin2, kout2, cube.shape[0])

    # find average continuum
    cont1 = np.nanmean(cube[kout1:kin1, :, :], axis=0)
    cont2 = np.nanmean(cube[kin2:kout2, :, :], axis=0)
    # Maybe don't use the continuum on one side (if contaminated)
    wt1, wt2 = usecont          
    avcont = (wt1*cont1 + wt2*cont2)/(wt1 + wt2)

    # Subtract to get pure line window
    linewin = cube[kin1:kin2, :, :] - avcont[None, :, :]
    # Now find moments
    nwin = kin2 - kin1
    _, _, winwavs = wcube.all_pix2world([0]*nwin, [0]*nwin, np.arange(kin1, kin2), 0)
    # Convert wavelengths to velocities
    winvels = waves2vels(winwavs*u.m.to(u.Angstrom), wav, helio_hdr, observatory='VLT4')
    mom0 = linewin.sum(axis=0)
    mom1 = np.sum(linewin*winvels[:, None, None], axis=0)
    vmean = mom1/mom0
    mom2 = np.sum(linewin*(winvels[:, None, None] - vmean)**2, axis=0)
    vsig = np.sqrt(mom2/mom0)
    maps = {'continuum': avcont, 'linesum': mom0,
            'mean': vmean.to(u.km/u.s).value, 'sigma': vsig.to(u.km/u.s).value}

    # Save an average spectrum as well
    if kout2 > cube.shape[0]:
        print('Out of bounds: truncating from', kout2, 'to', cube.shape[0])
        kout2 = cube.shape[0]
    nwinout = kout2 - kout1
    _, _, winwavsout = wcube.all_pix2world([0]*nwinout, [0]*nwinout, np.arange(kout1, kout2), 0)
    winvelsout = waves2vels(winwavsout*u.m.to(u.Angstrom), wav, helio_hdr, observatory='VLT4')

    # Just select region to SW of Trap
    xslice = slice(900, 1400)
    yslice = slice(400, 900)
    spec = {
        'vhel': winvelsout.to(u.km/u.s),
        'wav': winwavsout*u.m.to(u.Angstrom),
        'flux': np.nanmean(cube[kout1:kout2, yslice, xslice], axis=(1, 2)),
        'cont': np.nanmean(avcont[yslice, xslice])*np.ones(nwinout),
    }
    print(*[(k, len(v)) for k, v in spec.items()])
    return maps, spec
