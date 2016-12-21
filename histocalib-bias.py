from astropy.io import fits
from histocalib import histogram_calib_images, maxcount
MUSE_BIAS = 20.0
CONVERSION_FACTOR = 4.2e-6
if __name__ == '__main__':
    for f, vmax in maxcount.items():
        fn = 'muse-hr-cropspec1d-wfc3-{}.fits'.format(f)
        nv = fits.open(fn)[0].header['NAXIS1']
        muse_shift = MUSE_BIAS*CONVERSION_FACTOR*nv 
        print(histogram_calib_images(f, vmax, muse_shift))
