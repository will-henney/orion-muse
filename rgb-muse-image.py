TAB=[["C_II-723X-fluorescent", -3000, 31000, "linear"], ["LineMaps/linesum-Cl_IV-8046-multibin-SN0005", -800, 8200, "linear"], ["LineMaps/linesum-Ar_IV-4740-multibin-SN0005", 0, 9900, "linear"]]
SUFFIX="hi"
from astropy.io import fits
import aplpy

figfile = 'rgb-muse-{}.pdf'.format(SUFFIX)

# Unpack the channel info from the table
[rf, r1, r2, rs], [gf, g1, g2, gs], [bf, b1, b2, bs] = TAB
rgbfiles = [rf + '.fits', gf + '.fits', bf + '.fits']

# aplpy can only deal with the primary headers, so sort that out first
template = 'rgb-for-aplpy-{}.fits'
channels = ['red', 'green', 'blue']
newfiles = [template.format(chan) for chan in channels]
for newfn, fn in zip(newfiles, rgbfiles):
    hdu = fits.open(fn)['SCALED']
    hdu.writeto(newfn, clobber=True)

aplpy.make_rgb_image(newfiles, 'rgb-for-aplpy.png',
                     vmin_r=r1, vmin_g=g1, vmin_b=b1,
                     vmax_r=r2, vmax_g=g2, vmax_b=b2,
                     stretch_r=rs, stretch_g=gs, stretch_b=bs, 
                     make_nans_transparent=True)
f = aplpy.FITSFigure(newfiles[0])
f.show_rgb('rgb-for-aplpy.png')
# f.recenter(83.6875, -5.4167, width=0.25, height=0.167)
f.add_grid()
f.grid.set_color('white')
f.grid.set_alpha(0.2)
f.save(figfile)
