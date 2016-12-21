import sys
import os
import numpy as np
from astropy.io import fits
from matplotlib import pyplot as plt
import seaborn as sns

try: 
    infile = sys.argv[1]
except:
    sys.exit('Usage: {} FITSFILE'.format(sys.argv[0]))

try:
    YMIN, YMAX = float(sys.argv[2]), float(sys.argv[3])
except IndexError:
    YMIN, YMAX = 0.85, 1.15

basename = os.path.basename(infile)
baseroot, _ = os.path.splitext(basename)

hdu = fits.open(infile)[0]
if hdu.data is None:
    hdu = fits.open(infile)[1]
hdr = hdu.header

ny, nx = hdu.data.shape

# Size of chunks
mx, my = 290, 290
xchunks, ychunks = nx//mx, ny//my
xprofile = np.nanmean(hdu.data, axis=0)
yprofile = np.nanmean(hdu.data, axis=1)
plotfile = 'pattern-noise-xy-{}.pdf'.format(baseroot)
fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True)

x = np.arange(mx)
xpstack = []
for ichunk in range(xchunks):
    xp = xprofile[ichunk*mx:ichunk*mx + mx]
    ax1.plot(x, xp/xp.mean(), label=str(ichunk))
    xpstack.append(xp/xp.mean())
xpstack = np.vstack(xpstack)
xpm = np.median(xpstack, axis=0)
ax1.plot(x, xpm, 'k', lw=6, alpha=0.2)
ax1.legend(ncol=2)
ax1.set(
    ylim=[YMIN, YMAX],
    ylabel='x profile',
)

y = np.arange(my)
ypstack = []
for jchunk in range(ychunks):
    yp = yprofile[jchunk*my:jchunk*my + my]
    ax2.plot(y, yp/yp.mean(), label=str(jchunk))
    ypstack.append(yp/yp.mean())
ypstack = np.vstack(ypstack)
ypm = np.median(ypstack, axis=0)
ax2.plot(y, ypm, 'k', lw=6, alpha=0.2)
ax2.legend(ncol=2)
ax2.set(
    ylim=[YMIN, YMAX],
    xlabel='pixel in chunk',
    ylabel='y profile',
)

fig.set_size_inches(6, 6)
fig.tight_layout()
fig.savefig(plotfile)

# Finally, reconstruct the pattern noise image
patmap = np.ones_like(hdu.data)
tile = xpm[None, :]*ypm[:, None]
for jchunk in range(ychunks):
    yslice = slice(jchunk*my, jchunk*my + my)
    for ichunk in range(xchunks):
        xslice = slice(ichunk*mx, ichunk*mx + mx)
        patmap[yslice, xslice] = tile

fits.PrimaryHDU(data=patmap,
                header=hdr).writeto(infile.replace('.fits',
                                                   '-pattern.fits'), clobber=True)

fits.PrimaryHDU(data=hdu.data/patmap,
                header=hdr).writeto(infile.replace('.fits',
                                                   '-patfix.fits'), clobber=True)

print(plotfile, end='')
