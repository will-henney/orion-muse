from astropy.io import fits

filters_ = ["FQ575N", "FQ672N", "FQ674N", "F673N",
            "F487N", "F656N", "F658N", "F547M"]

def divide_fits_images(name1, name2, outname):
    hdu1 = fits.open(name1)[0]
    hdu2 = fits.open(name2)['DATA']
    fits.PrimaryHDU(header=hdu1.header, data=hdu1.data/hdu2.data).writeto(outname, clobber=True)
    print(outname)

if __name__ == '__main__':
    for f in filters_:
        divide_fits_images(
            'wfc3-new-resample-muse-{}.fits'.format(f),
            'muse-hr-new-cropimage-wfc3-{}.fits'.format(f),
            'wfc3-over-muse-new-calib-ratio-{}.fits'.format(f)
        )
