import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import aplpy
from astropy.io import fits
import pyregion
import numpy as np

Tsets = ['nominal', '2016-03-bias', '2016-03-allgrads']
ratios = ['5755-6583', '6716-6731', '4861-6563', '6716-6731-N']
vmin = {'5755-6583': 0.005, '6716-6731': 0.4, '4861-6563': 0.15}
vmax = {'5755-6583': 0.035, '6716-6731': 0.9, '4861-6563': 0.3}
cmap = matplotlib.cm.RdYlBu_r
ra0, dec0 = 83.80716, -5.403
dra, ddec = 0.025, 0.025
for r in ratios:
    rr = '-'.join(r.split('-')[:2])
    rlabel = rr.replace('-', ' / ')
    mfile = 'crop-ratio-{}-bin004.fits'.format(rr)
    for Tset in Tsets:
        wfiles = ['NebulioMUSE/wfc3-resample-{}-{}.fits'.format(Tset, r)
                  for Tset in Tsets]
        rfiles = ['NebulioMUSE/wfc3-over-muse-ratio-of-ratios-{}-{}.fits'.format(Tset, r)
                  for Tset in Tsets]
        for fitsfile in [mfile] + wfiles + rfiles:
            f = aplpy.FITSFigure(fitsfile) 
            f.recenter(ra0, dec0, width=dra, height=ddec)
            if 'ratio-of-ratios' in fitsfile:
                sweet_region = 'new-combined-sweet-black.reg'
                f.show_colorscale(interpolation='none', vmin=0.0, vmax=2.0, cmap=cmap)
            else:
                sweet_region = 'new-combined-sweet.reg'
                f.show_grayscale(interpolation='none', vmin=vmin[rr], vmax=vmax[rr])
            f.show_regions(sweet_region)
            f.add_label(0.2, 0.05, rlabel, 
                        horizontalalignment='left', relative=True)
            f.add_colorbar()
            figfile = fitsfile.replace('.fits', '.jpg')
            plt.gcf().set_size_inches(6, 5)
            plt.gcf().tight_layout()
            plt.gcf().savefig(figfile, dpi=300)
            print(figfile)
