import numpy as np
from astropy.io import fits
from pyneb_utils import rsii_T_den, rnii_T_den, rcliii_T_den, rsiii_T_den
from scipy.interpolate import griddata
import pyregion

# All values are now log10 until right at the end
Tmin, Tmax = 2.0, 5.0
Nmin, Nmax = 2.0, 5.0
nT, nN = 100, 100

Tgrid, Ngrid = np.linspace(Tmin, Tmax, nT), np.linspace(Nmin, Nmax, nN)
Nvalues, Tvalues = np.meshgrid(Ngrid, Tgrid)

rsii_grid = rsii_T_den(10**Tgrid, 10**Ngrid)
rnii_grid = rnii_T_den(10**Tgrid, 10**Ngrid)
points = np.array(zip(rsii_grid.ravel(), rnii_grid.ravel()))

rcliii_grid = rcliii_T_den(10**Tgrid, 10**Ngrid)
rsiii_grid = rsiii_T_den(10**Tgrid, 10**Ngrid)
points_iii = np.array(zip(rcliii_grid.ravel(), rsiii_grid.ravel()))

def T_den_from_rsii_rnii(rsii, rnii):
    """Uses grid interpolation to go from ratios to (Te, Ne)"""
    xi = np.array(zip(rsii, rnii))
    Te = griddata(points, Tvalues.ravel(), xi)
    Ne = griddata(points, Nvalues.ravel(), xi)
    # Repeat but with nearest neighbour
    Te_n = griddata(points, Tvalues.ravel(), xi, method='nearest')
    Ne_n = griddata(points, Nvalues.ravel(), xi, method='nearest')
    # And fill in any NaN or infinite values using these
    m = ~np.isfinite(Te)
    Te[m] = Te_n[m]
    m = ~np.isfinite(Ne)
    Ne[m] = Ne_n[m]
    return 10**Te, 10**Ne

def T_den_from_rcliii_rsiii(rcliii, rsiii):
    """Uses grid interpolation to go from ratios to (Te, Ne)"""
    xi = np.array(zip(rcliii, rsiii))
    Te = griddata(points_iii, Tvalues.ravel(), xi)
    Ne = griddata(points_iii, Nvalues.ravel(), xi)
    # Repeat but with nearest neighbour
    Te_n = griddata(points_iii, Tvalues.ravel(), xi, method='nearest')
    Ne_n = griddata(points_iii, Nvalues.ravel(), xi, method='nearest')
    # And fill in any NaN or infinite values using these
    m = ~np.isfinite(Te)
    Te[m] = Te_n[m]
    m = ~np.isfinite(Ne)
    Ne[m] = Ne_n[m]
    return 10**Te, 10**Ne


if __name__ == '__main__':

    hduA = fits.open("newratio-6716-6731-F673N.fits")[0]
    hduB = fits.open("newratio-5755-6583-F547M-deredden-2874.fits")[0]

    include = pyregion.open("will-nii-sweet-spot.reg")
    exclude = pyregion.open("will-nii-exclude.reg")
    m = include.get_mask(hdu=hduA) & (~exclude.get_mask(hdu=hduA))
    include = pyregion.open("will-sii-sweet-spot.reg")
    exclude = pyregion.open("will-sii-exclude.reg")
    m = m & include.get_mask(hdu=hduA) & (~exclude.get_mask(hdu=hduA))
    m = m & np.isfinite(hduA.data) & np.isfinite(hduB.data)
    Te = np.empty_like(hduA.data)
    Ne = np.empty_like(hduA.data)
    Te[m], Ne[m] = T_den_from_rsii_rnii(hduA.data[m], hduB.data[m])
    Te[~m], Ne[~m] = np.nan, np.nan

    fits.PrimaryHDU(header=hduA.header, data=Te).writeto('new-derived-Te.fits', clobber=True)
    fits.PrimaryHDU(header=hduA.header, data=Ne).writeto('new-derived-Ne.fits', clobber=True)
