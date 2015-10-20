from __future__ import print_function
import numpy as np
from astropy.io import fits

def reassemble(nsec=8, suffix='rebin05x05'):
    template = 'muse-hr-data-wavsec{isec}-{suffix}.fits'

    cubes = []
    for isec in range(nsec):
        fn = template.format(isec=isec, suffix=suffix)
        hdulist = fits.open(fn)
        hdu = hdulist[0]
        # The full cube will use the header from the first section
        if isec == 0:
            hdr = hdu.header.copy()
        # Save this section of the cube 
        cubes.append(hdu.data.astype(np.float32))
        print('Section', isec, 'saved')
        hdulist.close()
    newhdu = fits.PrimaryHDU(header=hdr, data=np.concatenate(cubes))
    outfile = 'muse-hr-fullcube-{suffix}.fits'.format(suffix=suffix)
    newhdu.writeto(outfile, clobber=True)
    return outfile

if __name__ == '__main__':

    print('Reassembling cube ...')
    print('... done:', reassemble())
