from astropy.io import fits

def save_ew(fname, hdu_id=0):
    '''Calculate and save EW for line map in file named `fname`

It is assumed that the continuum can be found by transforming the filename
'''
    hdu = fits.open(fname)[hdu_id]
    cname = fname.replace('linesum', 'continuum')
    cdata = fits.open(cname)[hdu_id].data
    hdu.data /= cdata
    ename = fname.replace('linesum', 'ew')
    hdu.writeto(ename, clobber=True)
