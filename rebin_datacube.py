from __future__ import print_function
import sys
import numpy as np
from astropy.io import fits

def rebin_xyimage(im, mx=5, my=5):
    ny, nx = im.shape
    # Shape of new rebinned array
    nny, nnx = ny//my, nx//mx
    # Shave a bit off original array so shape is multiple of m
    # And then concertina each axis to be 4-dimensional
    im4d = im[:nny*my, :nnx*mx].reshape((nny, my, nnx, mx))
    # Average along the mx, my axes
    return np.nansum(im4d, axis=(1, 3))


def rebin_hdu(hdu, m=(5, 5)):
    mx, my = m                  # FITS axis order
    nv, ny, nx = hdu.data.shape  # Python axis order
    nny, nnx = ny//my, nx//mx
    newdata = np.empty((nv, nny, nnx))
    for k in range(nv):
        print('Rebinning plane {}/{}'.format(k+1, nv))
        newdata[k] = rebin_xyimage(hdu.data[k], mx, my)
    newhdu = fits.PrimaryHDU(header=hdu.header, data=newdata)
    # New pixel deltas are bigger
    if 'CDELT1' in newhdu.header:
        newhdu.header['CDELT1'] *= mx
        newhdu.header['CDELT2'] *= my
    else:
        # assume no rotation
        newhdu.header['CD1_1'] *= mx
        newhdu.header['CD2_2'] *= my
      
         # pix=0.5 is left edge of first pixel
    newhdu.header['CRPIX1'] = 0.5 + (newhdu.header['CRPIX1'] - 0.5)/mx
    newhdu.header['CRPIX2'] = 0.5 + (newhdu.header['CRPIX2'] - 0.5)/my

    return newhdu

if __name__ == '__main__':
    try:
        infilename = sys.argv[1]
        m = int(sys.argv[2]), int(sys.argv[3])
    except IndexError:
        sys.exit('Usage: {} FITSFILE BINX BINY'.format(sys.argv[0]))

    hdu = fits.open(infilename)['DATA']
    newsuffix = '-rebin{:02d}x{:02d}.fits'.format(*m)
    outfilename = infilename.replace('.fits', newsuffix)
    rebin_hdu(hdu, m).writeto(outfilename, clobber=True)
