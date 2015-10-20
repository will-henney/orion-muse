import sys
import numpy as np
from astropy.io import fits
from astropy.wcs import WCS, WCSSUB_SPECTRAL
import astropy.units as u

def crop_muse_to_wfc3(fid):
    """Cut out a section of the MUSE image to match the WFC3 field"""
    wname = 'wfc3-resample-muse-{}.fits'.format(fid)
    mname = 'muse-hr-image-wfc3-{}.fits'.format(fid)
    whdu = fits.open(wname)[0]
    mhdu = fits.open(mname)['DATA']
    # Also get the spectral data cube multiplied by filter throughput
    shdu = fits.open(mname.replace('-image-', '-transwin-'))['DATA']
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

    # Finally, as a bonus, calculate the 1-D average spectrum from the cube
    spec = np.nanmean(shdu.data[:, start[1]:stop[1], start[0]:stop[0]], axis=(-1, -2))

    # Convert from 1e-20 flux-per-pixel to surface brightness units (flux per sr)
    pixel_area_sr = np.product(np.abs(wcs_m.wcs.cdelt))*(u.deg.to(u.radian))**2
    spec *= 1e-20/pixel_area_sr
    # extract only the spectral part of the cube's WCS
    wcs_s = WCS(shdu.header).sub([WCSSUB_SPECTRAL])
    oshdu = fits.PrimaryHDU(header=wcs_s.to_header(), data=spec)
    oshdu.header['BUNIT'] = 'erg/s/cm**2/sr/Angstrom'
    # Fix up the wavelngth units to angstrom
    oshdu.header['CDELT1'] *= 1e10
    oshdu.header['CRVAL1'] *= 1e10
    oshdu.header['CUNIT1'] = 'Angstrom'
    # And record the window from the MUSE full field that was extracted
    oshdu.header['MUSE_X1'] = start[0] + 1, 'Extracted window: start X pixel' 
    oshdu.header['MUSE_X2'] = stop[0] + 1, 'Extracted window: stop X pixel' 
    oshdu.header['MUSE_Y1'] = start[1] + 1, 'Extracted window: start Y pixel' 
    oshdu.header['MUSE_Y2'] = stop[0] + 1, 'Extracted window: stop Y pixel' 
    oshdu.writeto(mname.replace('-image-', '-cropspec1d-'), clobber=True)

    return oname


if __name__ == '__main__':
    try:
        filter_id = sys.argv[1]
    except:
        print('Usage:', sys.argv[0], 'FILTER')

    print(crop_muse_to_wfc3(filter_id))
