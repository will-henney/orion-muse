from astropy.io import fits
hduA = fits.open('LineMaps/continuum-N_II-5755-bin016.fits')['SCALED']
hduB = fits.open('LineMaps/continuum-Fe_III-5270-bin016.fits')['SCALED']
hduA.data /= hduB.data
hduA.writeto('cont-ratio-5755-to-5270-bin016.fits', clobber=True)
