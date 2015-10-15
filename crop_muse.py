import sys
import numpy as np
from astropy.io import fits
from astropy.wcs import WCS

def crop_muse_to_wfc3(fid):
    """Cut out a section of the MUSE image to match the WFC3 field"""
    wname = 'wfc3-resample-muse-{}.fits'.format(fid)
    mname = 'muse-hr-image-wfc3-{}.fits'.format(fid)
    whdu = fits.open(wname)[0]
    mhdu = fits.open(mname)['DATA']
    wcs_w = WCS(whdu.header).celestial
    wcs_m = WCS(mhdu.header).celestial
    # Check that the two images have the same reference point in RA, DEC
    assert np.all(wcs_w.wcs.crval == wcs_m.wcs.crval)
    # And that the pixel scales are the same
    assert np.all(wcs_w.wcs.cdelt == wcs_m.wcs.cdelt)
    assert np.all(wcs_w.wcs.pc == wcs_m.wcs.pc)

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
    start[0] += 1
    stop[0] += 1

    # Check that these are within bounds of the original MUSE grid
    assert np.all(start >= 0.0)
    assert np.all(stop < shape_m)

    # Crop the MUSE data array to the start:stop indices, remembering
    # that python axis order is backwards with respect to FITS axis
    # order
    mhdu.data = mhdu.data[start[1]:stop[1], start[0]:stop[0]]

    # And copy the WFC3 wcs into the new MUSE header
    mhdu.header.update(wcs_w.to_header())

    # Write out the new cropped MUSE image
    oname = mname.replace('-image-', '-cropimage-')
    mhdu.writeto(oname, clobber=True)

    return oname


if __name__ == '__main__':
    try:
        filter_id = sys.argv[1]
    except:
        print('Usage:', sys.argv[0], 'FILTER')

    print(crop_muse_to_wfc3(filter_id))
