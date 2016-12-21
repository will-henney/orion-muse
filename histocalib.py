from __future__ import print_function
import numpy as np
from astropy.io import fits
from astropy.convolution import convolve, Gaussian2DKernel
from matplotlib import pyplot as plt
import pyregion
import seaborn as sns
from specplot1d_utils import plot_1d_spec_from_fits

maxcount = {
    "fq575n": 0.4,
    "fq672n": 0.6,
    "fq674n": 0.75,
    "f673n" : 2.5,
    # "f469n" : 0.5,
    "f487n" : 7.0,
    "f656n" : 32.0, 
    "f658n" : 11.0, 
    "f547m" : 7.0, 
    # "f502n" : 20.0,
}
GAMMA = 2.0
T2DIR = '/Users/will/Work/RubinWFC3/Tsquared/'

# CORRECT_LINEAR_GRADIENT = ['fq672n', ]
CORRECT_LINEAR_GRADIENT = ["fq575n", "fq672n", "fq674n", "f673n",
                           "f487n" , "f656n" , "f658n" , "f547m"]
cmap = sns.light_palette((260, 50, 30), input="husl", as_cmap=True)
# cmap = plt.cm.gray_r

def histogram_calib_images(f, vmax=1.0, muse_shift=0.0):
    name1 = 'wfc3-new-resample-muse-{}.fits'.format(f)
    name2 = 'muse-hr-new-cropimage-wfc3-{}.fits'.format(f)
    if muse_shift == 0.0:
        pltname = 'wfc3-vs-muse-new-calib-{}.pdf'.format(f)
    else:
        pltname = 'wfc3-vs-muse-new-calib-{}-bias-shift.pdf'.format(f)

    hdu1 = fits.open(name1)[0]
    name3 = 'wfc3-over-muse-surface-calib-ratio-{}.fits'.format(f)
    if f in CORRECT_LINEAR_GRADIENT:
        # Linear correction surface to fix flat field
        hdus = fits.open(name3)[0]
        hdu1.data /= hdus.data
    hdu2 = fits.open(name2)['DATA']
    hduc = fits.open('wfc3-new-resample-muse-f547m.fits')[0]
    sweet_region = pyregion.open(T2DIR + 'sweet-3regions.reg')
    x, y = hdu2.data - muse_shift, hdu1.data
    xmin, xmax = ymin, ymax = 0.0, vmax
    ew = y/hduc.data
    # mask out silly values
    m = np.isfinite(x) & np.isfinite(y/x) & (np.abs(np.log10(y/x)) < 1.0)
    m = m & sweet_region.get_mask(hdu=hdu1)
    H, xedges, yedges = np.histogram2d(x[m], y[m], 100,
                                       [[xmin, xmax], [ymin, ymax]],
                                       weights=y[m])

    # Just take the median and weighted means of the ratio
    mm = m & (np.abs(np.log10(y/x)) < 0.3) # More restrictive mask
    median_slope = np.median(y[mm]/x[mm])
    weighted_mean_slope = np.average(y[mm]/x[mm], weights=y[mm])
    # And then take mean +/- std of those two to give the final answer
    slope = np.mean([median_slope, weighted_mean_slope])
    dslope = abs(median_slope - weighted_mean_slope)
    # Convert to a polynomial: linear with zero intercept
    pbest = np.poly1d([slope, 0.0])

    # Distribution of ratio, divided into two subsamples on EW or counts
    ratio = y/x
    # if f == 'f547m':
    if True:
        # Divide into high and low continuum counts
        s = 'counts'
        mlo = (y < y[mm].mean()) & mm
        mhi = (y >= y[mm].mean()) & mm
    else:
        # Divide into high and low EW
        s = f.upper() + '/F547M'
        mlo = (ew < np.median(ew[mm])) & mm
        mhi = (ew >= np.median(ew[mm])) & mm
    assert mlo.sum() > 0, f

    # Separate averages for the two subsamples
    ratio_hi = np.average(ratio[mhi], weights=y[mhi])
    ratio_lo = np.average(ratio[mlo], weights=y[mlo])
    # Difference in subsamples ... 
    ddslope = 0.5*np.abs(ratio_hi - ratio_lo)
    # ... which feeds in to the uncertainty in slope
    dslope = max(dslope, ddslope)

    # H = convolve(H, Gaussian2DKernel(1.0))
    fitcolor = (1.0, 0.5, 0.0)
    fitlabel = "y = ({:.3f} +/- {:.3f}) x".format(slope, dslope)
    fig, ax = plt.subplots(1, 1)
    ax.imshow((H.T)**(1.0/GAMMA), extent=[xmin, xmax, ymin, ymax],
              interpolation='none', aspect='auto', origin='lower', 
              cmap=cmap, alpha=1.0)
    ax.plot([0.0, x[m].max()], [0.0, x[m].max()], '-', alpha=1.0,
            lw=1, c='w', label=None)
    ax.plot([0.0, x[m].max()], pbest([0.0, x[m].max()]), '-', alpha=1.0,
            lw=1, c=fitcolor, label=fitlabel)

    leg = ax.legend(loc='upper left', title='Best Fit', frameon=True, fancybox=True)
    leg.get_title().set_fontsize('small')
    ax.set_ylabel(
        'Observed WFC3 {} count rate, electron/s/pixel'.format(f.upper()))
    ax.set_xlabel(
        'MUSE-predicted WFC3 {} count rate, electron/s/pixel'.format(f.upper()))
    ax.set_xlim(xmin, xmax)
    ax.set_ylim(ymin, ymax)


    # Now do 1-D histogram of the ratios 
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
    ax2.set_xlabel('Observed / Predicted', fontsize='xx-small')
    ax2.set_ylabel('Weighted PDF Histograms', fontsize='xx-small')
#    ax2.set_title('PDF', fontsize='x-small')

# inset axis at the bottom right
    ax3 = fig.add_axes([0.6, 0.2, 0.3, 0.3])
    fn = 'muse-hr-cropspec1d-wfc3-{}.fits'.format(f)
    plot_1d_spec_from_fits(fn, ax3, fontsize='xx-small')
    ax3.tick_params(labelsize='xx-small')
    ax3.set_title(f.upper())


    fig.set_size_inches(6, 6)
    fig.savefig(pltname)

    return [pltname, fitlabel]


if __name__ == '__main__':
    for f, vmax in maxcount.items():
        print(histogram_calib_images(f, vmax))
