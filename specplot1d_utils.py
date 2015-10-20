from __future__ import print_function
import sys
from astropy.io import fits
from astropy.wcs import WCS
from astropy import units as u
from matplotlib import pyplot as plt
import seaborn as sns

def plot_1d_spec_from_fits(fn, ax, fontsize=None):
    """Plots spectrum from filename `fn` onto pre-existing axis `ax`"""
    hdu = fits.open(fn)[0]
    spec = hdu.data/1e-3
    w = WCS(hdu.header)
    nwav = len(spec)
    wavs, = w.all_pix2world(range(nwav), 0)
    wavs *= u.m.to(u.Angstrom)
    #ax.plot(wavs, spec, drawstyle='steps-mid')
    ax.bar(wavs, spec, align='center', linewidth=0)
    ax.set_xlim(wavs.min(), wavs.max())
    ax.set_xlabel('Observed Air Wavelength, Angstrom', fontsize=fontsize)
    ax.set_ylabel('Filter Throughput x Brightness\n 0.001 erg/s/cm^2/sr/Angstrom', fontsize=fontsize)


if __name__ == '__main__':
    try:
        filt = sys.argv[1]
    except IndexError:
        print('Usage:', sys.argv[0], 'FILTER')
    fig, ax = plt.subplots(1, 1)
    fn = 'muse-hr-cropspec1d-wfc3-{}.fits'.format(filt)
    plot_1d_spec_from_fits(fn, ax)
    # ax.set_yscale('log')
    # ax.set_ylim(1e-7, None)
    fig.savefig(sys.argv[0].replace('.py', '-test-{}.pdf'.format(filt)))
