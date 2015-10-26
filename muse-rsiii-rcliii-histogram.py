from __future__ import print_function
import sys
import numpy as np
from astropy.io import fits
from matplotlib import pyplot as plt
import pyregion
t2dir = '/Users/will/Work/RubinWFC3/Tsquared'
sys.path.append(t2dir)
from pyneb_utils import rcliii_T_den, rsiii_T_den

try:
    region = sys.argv[1]
except IndexError:
    sys.exit('Usage: {} REGION [SUFFIX]'.format(sys.argv[0]))

try:
    suffix = '-' + sys.argv[2]
except IndexError:
    suffix = ''

ssiii_hdu = fits.open("Linemaps/linesum-S_III-9069{}.fits".format(suffix))['SCALED']
hduA = fits.open("Linemaps/ratio-5538-5518{}.fits".format(suffix))['SCALED']
hduB = fits.open("Linemaps/ratio-6312-9069-deredden{}.fits".format(suffix))['SCALED']
# hduB = fits.open(prefix + "ratio-5755-6583.fits")[0]

xmin, xmax, ymin, ymax = 0.5, 2.5, 0.02, 0.08

m = np.isfinite(hduA.data) & np.isfinite(hduB.data)
if 'sweet' in region.lower():
    title = 'Orion S'
    include = pyregion.open(t2dir + '/will-nii-sweet-spot-wcs.reg')
    exclude = pyregion.open(t2dir + '/will-nii-exclude-wcs.reg')
    m = m & include.get_mask(hdu=ssiii_hdu) & (~exclude.get_mask(hdu=ssiii_hdu))
     
    include = pyregion.open(t2dir + '/will-sii-sweet-spot-wcs.reg')
    exclude = pyregion.open(t2dir + '/will-sii-exclude-wcs.reg')
    m = m & include.get_mask(hdu=ssiii_hdu) & (~exclude.get_mask(hdu=ssiii_hdu))
else:
    title = 'Full nebula'

x = hduA.data[m]
y = hduB.data[m]
w = ssiii_hdu.data[m]
gamma = 1.5
fig, ax = plt.subplots(1, 1)
H, xedges, yedges = np.histogram2d(x, y, 50,
                                   [[xmin, xmax], [ymin, ymax]],
                                   weights=w)
ax.imshow((H.T)**(1.0/gamma), extent=[xmin, xmax, ymin, ymax],
          interpolation='nearest', aspect='auto', origin='lower', 
          cmap=plt.cm.gray_r, alpha=1.0)

denrange = np.linspace(0.0, 8.e4, 200)
for tem in [7.5e3, 8e3, 8.5e3, 9e3, 9.5e3, 1e4]:
    ax.plot(rcliii_T_den(tem, denrange), rsiii_T_den(tem, denrange),
            label="T = {:.0f} K".format(tem))

temrange = np.linspace(0.0, 20000.0, 100)
for den in [1000, 2000, 4000, 8000, 16000, 32000]:
    ax.plot(rcliii_T_den(temrange, den), rsiii_T_den(temrange, den),
            '--', label="n = {:.0f} pcc".format(den))

ax.set_xlim(xmin, xmax)
ax.set_ylim(ymin, ymax)
ax.set_xlabel('[Cl III] 5538 / 5518')
ax.set_ylabel('Dereddened [S III] 6300 / 9069')
ax.grid()
ax.legend(ncol=2, fontsize='x-small', handlelength=2.2, numpoints=1)

fig.set_size_inches(7, 7)
pltfile = 'muse-rsiii-rcliii-histogram-{}{}.pdf'.format(region, suffix)


fig.savefig(pltfile)
print(pltfile)
