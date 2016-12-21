from astropy.io import fits
hduA = fits.open('LineMaps/linesum-C_II-7236-multibin-SN0010.fits')['SCALED']
hduAA = fits.open('LineMaps/linesum-C_II-7231-multibin-SN0010.fits')['SCALED']
hduB = fits.open('LineMaps/linesum-He_I-6678-multibin-SN0005.fits')['SCALED']
hduA.data += hduAA.data
hduA.data -= hduB.data/13.0
hduA.writeto('C_II-723X-fluorescent.fits', clobber=True)
