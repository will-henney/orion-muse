from astropy.io import fits
hduA = fits.open('LineMaps/continuum-N_II-5755-bin004.fits')['SCALED']
hduB = fits.open('muse-hr-image-wfc3-f547m-bin004.fits')['SCALED']
hduA.data /= hduB.data
hduA.writeto('cont-ratio-5755-to-F547M-bin004.fits', clobber=True)
