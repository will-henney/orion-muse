from astropy.io import fits
datadir = "/Volumes/SSD-1TB/OrionMuse/LineMaps"
hduA = fits.open(f'{datadir}/linesum-C_II-7236-bin001.fits')['SCALED']
hduAA = fits.open(f'{datadir}/linesum-C_II-7231-bin001.fits')['SCALED']
hduB = fits.open(f'{datadir}/linesum-He_I-6678-bin001.fits')['SCALED']
hduA.data += hduAA.data
hduA.data -= hduB.data/13.0
hduA.writeto('C_II-723X-fluorescent-bin001.fits', clobber=True)
