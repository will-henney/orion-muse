from fig_utils import fig_ax_im_from_fits


quantities = {
    'muse-derived-Te': (6500, 14000),
    'muse-derived-Te-iii': (6500, 14000),
    'muse-derived-Ne': (0.0, 2.0e4),
    'muse-derived-Ne-iii': (0.0, 2.0e4),
}

for SN in 3, 10, 30, 100, 300:
    for quant, (vmin, vmax) in quantities.items():
        fn = ('LineMaps/{}-multibin-SN{:04d}.fits'.format(quant, SN))
        try:
            fig, ax, im = fig_ax_im_from_fits(fn, vmin=vmin, vmax=vmax, cmap=plt.cm.gray_r)
        except IOError:
            continue
        figfile = fn.replace('.fits', '.jpg')
        fig.set_size_inches(8, 6)
        fig.savefig(figfile, dpi=200)
