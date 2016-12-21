from astropy.io import fits

ratios = ['5755-6583', '6716-6731', '4861-6563', '6716-6731-N']

def divide_fits_images(name1, name2, outname):
    hdu1 = fits.open(name1)[0]
    hdu2 = fits.open(name2)['SCALED']
    fits.PrimaryHDU(header=hdu1.header,
                    data=hdu1.data/hdu2.data).writeto(outname, clobber=True)
    print(outname)

if __name__ == '__main__':
    for Tset in 'nominal', '2016-03-bias', '2016-03-allgrads':
        for r in ratios:
            # remove third element for MUSE ratio file
            rr = '-'.join(r.split('-')[:2]) 
            divide_fits_images(
                'NebulioMUSE/wfc3-resample-{}-{}.fits'.format(Tset, r),
                'crop-ratio-{}-bin004.fits'.format(rr),
                'NebulioMUSE/wfc3-over-muse-ratio-of-ratios-{}-{}.fits'.format(Tset, r)
        )
