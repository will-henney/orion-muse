from __future__ import print_function
import sys
from astropy.io import fits
from astropy import wcs
from astropy.table import Table
import pysynphot
import numpy as np
from astropy import units as u
def air_refractive_index(wav):
    """Equation (65) of Greisen et al 2006 for the refractive index of air
at STP.  Input wavelength 'wav' should be in microns or in any
'astropy.units' unit. It does not matter if 'wav' is on air or vacuum
scale

    """
    try:
        # Convert to microns if necessary
        wavm = wav.to(u.micron).value
    except AttributeError:
        # Assume already in microns
        wavm = wav
    return 1.0 + 1e-6*(287.6155 + 1.62887/wavm**2 + 0.01360/wavm**4)

def bp_fullname(instrument, filter_):
    if instrument.lower() == 'wfc3':
        return 'wfc3,uvis1,'+filter_.lower()
    elif instrument.lower() == 'acs':
        return 'acs,wfc1,'+filter_.lower()
    elif instrument.lower() == 'wfpc2':
        return 'wfpc2,'+filter_.lower()
    else:
        raise NotImplementedError('Unknown instrument: ' + instrument)
from astropy import units as u
WFC3_CONSTANT = 0.0840241
MUSE_FLUX_UNITS = 1e-20 
MUSE_PIXEL_AREA_SR = (0.2*u.arcsec).to(u.radian)**2


def bandpass_flatten(instrument, bpname):
    filename = 'muse-hr-window-{}-{}.fits'.format(instrument, bpname)
    hdulist = fits.open(filename)
    hdu = hdulist['DATA']
    w = wcs.WCS(hdu.header)
    NV, NY, NX = hdu.data.shape
    # construct array of observed air wavelengths (at image center to be safe)
    _, _, wavs = w.all_pix2world([NX/2]*NV, [NY/2]*NV, np.arange(NV), 0) 
    # Make dimensional
    wavs *= u.m
    # Convert to vacuum scale
    wavs *= air_refractive_index(wavs)

    # Get bandpass for filter
    fn = bp_fullname(instrument, bpname)
    bp = pysynphot.ObsBandpass(fn)
    # Calculate transmission curve at the observed wavelengths
    T = bp(wavs.to(u.Angstrom).value)
    # Weight by transmission curve and save that
    hdu.data *= T[:, None, None]
    hdulist.writeto(filename.replace('-window-', '-transwin-'), clobber=True)
    # Integrate over wavelength, already weighted by transmission curve. But
    # now need to multiply by wavelength, put in brightness units, and
    # convert to WFC3 electron/s/pixel
    hdu.data *= WFC3_CONSTANT*MUSE_FLUX_UNITS/MUSE_PIXEL_AREA_SR
    hdu.data *= wavs.to(u.Angstrom).value[:, None, None]
    hdu.data = hdu.header['CDELT3']*np.sum(hdu.data, axis=0)
    hdu.header['BUNIT'] = 'electron/s/(0.03962 arcsec)**2'
    hdulist.writeto(filename.replace('-window-', '-image-'), clobber=True)

if __name__ == '__main__':
    try:
        instrument, bpname = sys.argv[1:]
        bandpass_flatten(instrument, bpname)
    except IndexError:
        print('Usage:', sys.argv[0], 'INSTRUMENT FILTER')
