from __future__ import print_function
import numpy as np
from astropy.io import fits
from astropy.convolution import convolve, Gaussian2DKernel
from matplotlib import pyplot as plt
import seaborn as sns
from specplot1d_utils import plot_1d_spec_from_fits

maxcount = {
    "fq575n": 0.4,
    "fq672n": 0.6,
    "fq674n": 0.75,
    "f673n" : 2.5,
    "f469n" : 0.5,
    "f487n" : 10.0,
    "f656n" : 40.0, 
    "f658n" : 11.0, 
    "f547m" : 7.0, 
    "f502n" : 20.0,
}
GAMMA = 2.0

cmap = sns.light_palette((260, 50, 30), input="husl", as_cmap=True)
# cmap = plt.cm.gray_r

def histogram_calib_images(f, vmax=1.0):
    name1 = 'wfc3-resample-muse-{}.fits'.format(f)
    name2 = 'muse-hr-cropimage-wfc3-{}.fits'.format(f)
    pltname = 'wfc3-vs-muse-calib-{}.pdf'.format(f)
    hdu1 = fits.open(name1)[0]
    hdu2 = fits.open(name2)['DATA']
    hduc = fits.open('wfc3-resample-muse-f547m.fits')[0]
    x, y = hdu2.data, hdu1.data
    xmin, xmax = ymin, ymax = 0.0, vmax
    ew = y/hduc.data
    # mask out silly values
    m = np.isfinite(x) & np.isfinite(y/x) & (np.abs(np.log10(y/x)) < 1.0)
    H, xedges, yedges = np.histogram2d(x[m], y[m], 200,
                                       [[xmin, xmax], [ymin, ymax]],
                                       weights=y[m])
    # Fit a straight line
    mm = m & (x > 0.05*xmax) & (y > 0.05*ymax) & (x < 0.5*xmax) & (y < 0.5*ymax) & (np.abs(np.log10(y/x)) < 0.3)
    # First, linear fit y(x) = m x + c
    y_x_linfit = np.polyfit(x[mm], y[mm], 1, w=y[mm])
    # Second, linear fit x(y) = m y + c
    x_y_linfit = np.polyfit(y[mm], x[mm], 1, w=y[mm])
    # Convert this from x(y) -> y(x)
    # If x = m y + c, then y = (1/m) x - c/m
    y_x_altfit = np.array([1./x_y_linfit[0], -x_y_linfit[1]/x_y_linfit[0]])
    # Take average and spread of these two fits
    y_x_bestfit = 0.5*(y_x_linfit + y_x_altfit)
    y_x_errfit = 0.5*np.abs(y_x_linfit - y_x_altfit)

    pbest =  np.poly1d(y_x_bestfit)

    # H = convolve(H, Gaussian2DKernel(1.0))
    fitcolor = (1.0, 0.5, 0.0)
    fitlabel = "y = ({:.2f} +/- {:.2f}) x + ({:.2f} +/- {:.2f})".format(
        y_x_bestfit[0], y_x_errfit[0], y_x_bestfit[1], y_x_errfit[1])
    fig, ax = plt.subplots(1, 1)
    ax.imshow((H.T)**(1.0/GAMMA), extent=[xmin, xmax, ymin, ymax],
              interpolation='none', aspect='auto', origin='lower', 
              cmap=cmap, alpha=1.0)
    ax.plot([0.0, x[m].max()], [0.0, x[m].max()], '-', alpha=1.0,
            lw=1, c='w', label=None)
    ax.plot([0.0, x[m].max()], pbest([0.0, x[m].max()]), '-', alpha=1.0,
            lw=1, c=fitcolor, label=fitlabel)

    leg = ax.legend(loc='upper left', title='Linear Fit', frameon=True, fancybox=True)
    leg.get_title().set_fontsize('small')
    ax.set_ylabel(
        'Observed WFC3 {} count rate, electron/s/pixel'.format(f.upper()))
    ax.set_xlabel(
        'MUSE-predicted WFC3 {} count rate, electron/s/pixel'.format(f.upper()))
    ax.set_xlim(xmin, xmax)
    ax.set_ylim(ymin, ymax)

    # Now do 1-D histogram of the deviations from the model
    ratio = y/pbest(x)
    if f == 'f547m':
        # Divide into high and low continuum counts
        s = 'counts'
        mlo = (y < y[m].mean()) & m
        mhi = (y >= y[m].mean()) & m
    else:
        # Divide into high and low EW
        s = f.upper() + '/F547M'
        mlo = (ew < np.median(ew[m])) & m
        mhi = (ew >= np.median(ew[m])) & m

    assert mlo.sum() > 0, f
    
    # inset axis at the top left
    ax2 = fig.add_axes([0.2, 0.55, 0.25, 0.25])
    ax2.hist(ratio[mlo], bins=100, range=(0.5, 1.5),
             normed=True, weights=y[mlo], alpha=0.7, label='Low '+s)
    ax2.hist(ratio[mhi], bins=100, range=(0.5, 1.5),
             normed=True, weights=y[mhi], color='red', alpha=0.3, label='High '+s)
    ax2.set_xlim(0.5, 1.5)
    # leave more space at top
    y1, y2 = ax2.get_ylim()
    y2 *= 1.2
    ax2.set_ylim(y1, y2)
    ax2.legend(loc='upper left', fontsize='xx-small')
    ax2.tick_params(labelleft=False, labelsize='x-small')
    ax2.set_xlabel('(Observed Counts) / (Linear Fit)', fontsize='xx-small')
    ax2.set_ylabel('Weighted PDF Histograms', fontsize='xx-small')
#    ax2.set_title('PDF', fontsize='x-small')

    # inset axis at the bottom right
    ax3 = fig.add_axes([0.6, 0.2, 0.3, 0.3])
    fn = 'muse-hr-cropspec1d-wfc3-{}.fits'.format(f)
    plot_1d_spec_from_fits(fn, ax3, fontsize='xx-small')
    ax3.tick_params(labelsize='xx-small')
    ax3.set_title(f.upper())


    fig.set_size_inches(7, 7)
    fig.savefig(pltname)

    return [pltname, fitlabel]


if __name__ == '__main__':
    for f, vmax in maxcount.items():
        print(histogram_calib_images(f, vmax))
