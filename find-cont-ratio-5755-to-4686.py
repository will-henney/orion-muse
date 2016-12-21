from astropy.io import fits
hduA = fits.open('LineMaps/continuum-N_II-5755-bin004.fits')['SCALED']
hduB = fits.open('LineMaps/continuum-H_I-4686-bin004.fits')['SCALED']
hduA.data /= hduB.data
hduA.writeto('cont-ratio-5755-to-4686-bin004.fits', clobber=True)
