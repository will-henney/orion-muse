import sys
import numpy as np
from astropy.io import fits
from astropy.wcs import WCS, WCSSUB_SPECTRAL
import astropy.units as u

def crop_muse_to_wfc3(
        fid,
        muse_pattern='LineMaps/ratio-{}-bin004.fits',
        wfc3_pattern='NebulioMuse/wfc3-resample-2016-03-{}.fits',
        crop_pattern='crop-ratio-{}-bin004.fits',
        wfc3_hduname=0, muse_hduname='SCALED'
        ):
    """Cut out a section of the MUSE image to match the WFC3 field"""
    wname = wfc3_pattern.format(fid)
    mname = muse_pattern.format(fid)
    whdu = fits.open(wname)[wfc3_hduname]
    mhdu = fits.open(mname)[muse_hduname]
    wcs_w = WCS(whdu.header).celestial
    wcs_m = WCS(mhdu.header).celestial
    # Check that the two images have the same reference point in RA, DEC
    assert np.all(wcs_w.wcs.crval == wcs_m.wcs.crval)
    # And that the pixel scales are the same
    wcd = wcs_w.wcs.cdelt*wcs_w.wcs.pc
    mcd = wcs_m.wcs.cdelt*wcs_m.wcs.pc
    assert np.all(wcd == mcd), '{} differs from {}'.format(repr(wcd), repr(mcd)) 

    # The shapes of the two grids: (nx, ny) in FITS axis order
    shape_w = np.array([whdu.header['NAXIS1'], whdu.header['NAXIS2']])
    shape_m = np.array([mhdu.header['NAXIS1'], mhdu.header['NAXIS2']])

    # The difference in CRPIX values tells us the start indices (i, j)
    # for the crop window on the MUSE grid. Note that this is in
    # zero-based array indices
    start = wcs_m.wcs.crpix - wcs_w.wcs.crpix
    # The stop indices for the crop window 
    stop = start + shape_w

    # Shift 1 pixel to the right to do a coarse alignment correction
    start[0] += -1
    stop[0] += -1

    # Check that these are within bounds of the original MUSE grid
    assert np.all(start >= 0.0)
    assert np.all(stop < shape_m)

    # Crop the MUSE data array to the start:stop indices, remembering
    # that python axis order is backwards with respect to FITS axis
    # order
    newview = (slice(start[1], stop[1]), slice(start[0], stop[0]))
    mhdu.data = mhdu.data[newview]

    # And crop the WCS object in the same way
    wcs_new = wcs_m.slice(newview)
    # And copy the new cropped wcs into the new MUSE header
    mhdu.header.update(wcs_new.to_header())

    # Write out the new cropped MUSE image
    oname = crop_pattern.format(fid)
    mhdu.writeto(oname, clobber=True)

    return oname


if __name__ == '__main__':
    try:
        ratio_id = sys.argv[1]
    except:
        print('Usage:', sys.argv[0], 'RATIO')

    print(crop_muse_to_wfc3(ratio_id))
