import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import aplpy
from astropy.io import fits
import pyregion
import numpy as np

filters_ = ["FQ575N", "FQ672N", "FQ674N", "F673N",
            "F487N", "F656N", "F658N", "F547M"]
cmap = matplotlib.cm.RdYlBu_r
ra0, dec0 = 83.80716, -5.403
dra, ddec = 0.025, 0.025
sweet_region = 'new-combined-sweet.reg'
for filter_ in filters_:
    fitsfile = 'wfc3-over-muse-new-calib-ratio-{}.fits'.format(filter_) 
    f = aplpy.FITSFigure(fitsfile) 
    f.recenter(ra0, dec0, width=dra, height=ddec)
    f.show_colorscale(interpolation='none', vmin=0.5, vmax=1.5, cmap=cmap)
    for goodbad in 'good', 'bad':
        try:
            regfile = '{}-{}.reg'.format(filter_, goodbad)
            f.show_regions(regfile)
        except IOError:
            pass
    f.add_label(0.1, 0.1, filter_, relative=True)
    f.add_colorbar()
    f.colorbar.set_axis_label_text('WFC3 / MUSE')
    f.axis_labels.hide()
    f.tick_labels.hide()
    figfile = fitsfile.replace('.fits', '.jpg')
    plt.gcf().set_size_inches(6, 5)
    plt.gcf().savefig(figfile, dpi=300)
    print(figfile)

# Now plot all the regions on top of one another
fitsfile = 'wfc3-over-muse-new-calib-ratio-{}.fits'.format('F673N')
hdu = fits.open(fitsfile)[0]
hdu.data[np.isfinite(hdu.data)] = 1.0
hdu.data[~np.isfinite(hdu.data)] = 0.0
badregions = []
for filter_ in filters_:
    try:
        goodregion = pyregion.open(filter_ + '-good.reg')
        mask = goodregion.get_mask(hdu=hdu)
        hdu.data[mask] += 1.0
    except IOError:
        pass
for filter_ in filters_:
    try:
        badregion = pyregion.open(filter_ + '-bad.reg')
        mask = badregion.get_mask(hdu=hdu)
        hdu.data[mask] *= 0.3
        badregions.append(badregion)
    except IOError:
        pass
f = aplpy.FITSFigure(hdu)
f.recenter(ra0, dec0, width=dra, height=ddec)
f.show_grayscale(interpolation='none', invert=True)
f.add_colorbar()
#for region in badregions:
#    f.show_regions(region)
f.show_regions(sweet_region)
f.axis_labels.hide()
f.tick_labels.hide()
figfile = 'wfc3-over-muse-new-calib-ratio-SWEET.jpg'
plt.gcf().set_size_inches(6, 5)
plt.gcf().savefig(figfile, dpi=300)
print(figfile)
