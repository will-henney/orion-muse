from __future__ import print_function
import sys
import numpy as np
from astropy.table import Table
from astropy.io import fits
from astropy.wcs import WCS
from astropy import units as u
from misc_utils import sanitize_string
sys.path.append('/Users/will/Dropbox/OrionWest/')
from helio_utils import waves2vels

linetab = Table.read('basic-line-list.tab', format='ascii.tab')
wavsectab = Table.read('wavsec-startwavs.tab', format='ascii.tab')

def find_wavsec(wav):
    """Which wavsec chunk is this wavelength in"""
    return wavsectab['Section'][wav > wavsectab['CRVAL3']].max() 


def save_linemap_files(wav, species, mapdir='LineMaps', usecont=[1, 1]):
    wavsec = find_wavsec(wav)
    wavsec_blue = find_wavsec(wav - 8.0)
    wavsec_red = find_wavsec(wav + 8.0)
    if (wavsec_blue != wavsec) or (wavsec_red != wavsec):
        print('Uh, oh - line straddles wavsecs', wavsec_blue, wavsec, wavsec_red)
        if 5000.0 < wav < 6000.0:
            print('But luckily we can use F547M window instead')
            fn = 'muse-hr-window-wfc3-f547m.fits'
        else:
            sys.exit('Unable to find a cube that contains wavelength range - sorry!')
    else:
        fn = 'muse-hr-data-wavsec{}.fits'.format(wavsec)
    hdulist = fits.open(fn)
    hdu = hdulist['DATA']
    w = WCS(hdu.header)
    helio_hdr = fits.open('muse-hr-window-wfc3-f656n.fits')[0].header
    maps, spec =  extract_line_maps(wav, hdu.data, w, helio_hdr, usecont)
    wavid = str(int(wav+0.5))
    # Save the maps to FITS file 
    for mapid, mapdata in maps.items():
        mhdu = fits.PrimaryHDU(header=w.celestial.to_header(), data=mapdata)
        mname = '{}/{}-{}-{}.fits'.format(mapdir, mapid, species, wavid)
        mhdu.writeto(mname, clobber=True)
    sname = '{}/spec1d-{}-{}.tab'.format(mapdir, species, wavid)
    # And save the spectrum to TSV file
    Table(spec).write(sname, format='ascii.tab')

def extract_line_maps(wav, cube, wcube, helio_hdr,
                      usecont=[1, 1], dwav=4.0, dwavcont=8.0):
    """Extract line moments and continuum"""

    # outer range is wav +/- dwavcont
    wavs = (wav + dwavcont*np.array([-1, 1]))*u.Angstrom.to(u.m)
    _, _, (kout1, kout2) = wcube.all_world2pix([0, 0], [0, 0], wavs, 0)
    kout1, kout2 = int(kout1), int(kout2) + 2

    # inner range is wav +/- dwav
    wavs = (wav + dwav*np.array([-1, 1]))*u.Angstrom.to(u.m)
    _, _, (kin1, kin2) = wcube.all_world2pix([0, 0], [0, 0], wavs, 0)
    kin1, kin2 = int(kin1), int(kin2) + 2

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

    return maps, spec


for row in linetab:
    print(row['Ion'], row['wav0'])
    save_linemap_files(row['wav0'], sanitize_string(row['Ion']),
                       usecont=[row['blue cont'], row['red cont']])
