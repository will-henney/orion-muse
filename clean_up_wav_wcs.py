import sys
from astropy.io import fits
def clean_up_wav_wcs(filename):
    hdulist = fits.open(filename, mode='update')
    for hdu in hdulist:
        if hdu.header.get('CUNIT3') == 'm':
            # Change to Angstrom
            hdu.header['PC3_3'] *= 1e10
            hdu.header['CRVAL3'] *= 1e10
            hdu.header['CUNIT3'] = 'Angstrom'
            # And move scales to CDELT
            for i in '123':
                CDELTi = 'CDELT'+i
                # Sanity check
                assert hdu.header.get(CDELTi) == 1.0
                PCi_j = 'PC{0}_{0}'.format(i)
                hdu.header[CDELTi], hdu.header[PCi_j] = hdu.header[PCi_j], hdu.header[CDELTi] 
    hdulist.flush()


if __name__ == '__main__':
    try:
        fn = sys.argv[1]
        clean_up_wav_wcs(fn)
    except IndexError:
        print('Usage:', sys.argv[0], 'FITSFILE')
