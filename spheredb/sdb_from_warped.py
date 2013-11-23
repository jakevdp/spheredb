import numpy as np
from scipy.sparse import coo_matrix

from astropy.io import fits
from astropy.wcs import WCS
from astropy.units import Unit


def sdb_from_warped(fitsfile, hdunum=1):
    """Construct a sparse matrix representation of the FITS data

    Parameters
    ----------
    fitsfile : str
        A FITS file, which must be in HPX format.
    hdunum : int
        The number of the HDU to use (default = 1)

    Retrurns
    --------
    arr : scipy.sparse.coo_matrix
        The sparse representation of the data in HPX format
    """
    hdulist = fits.open(fitsfile)
    wcs = WCS(hdulist[hdunum].header)
    img = hdulist[hdunum].data

    # Check that unit is what we expect
    if wcs.wcs.ctype[0] != 'RA---HPX' or wcs.wcs.ctype[1] != 'DEC--HPX':
        raise ValueError("Expected HEALPix input. "
                         "Got {0}".format(wcs.wcs.ctype))

    if wcs.wcs.cunit[0] != Unit('deg') or wcs.wcs.cunit[1] != Unit('deg'):
        raise ValueError("Expected units of degrees.  "
                         "Got {0}".format(wcs.wcs.cunit))

    # Find cdelt along both axes
    try:
        cd = wcs.wcs.cd
    except:
        cd = None

    if cd is None:
        cdelt = w.cdelt
    elif cd.shape != (2, 2):
        raise ValueError('cd not understood')
    elif cd[0, 1] != 0 or cd[1, 0] != 0:
        raise ValueError('non-diagonal cd matrix not supported')
    else:
        cdelt = cd.diagonal()

    # Construct the matrix
    Nx = wcs.naxis1
    Ny = wcs.naxis2
    Nx_tot = int(np.round(360. / cdelt[0]))
    Ny_tot = int(np.round(180. / cdelt[1]))

    # sanity check for img shape
    assert img.shape == (Ny, Nx)

    # x counts starting at RA=0, y counts starting at DEC=-90
    x0, y0 = np.round(-wcs.wcs.crpix).astype(int)
    y0 += Ny_tot / 2

    ix, iy = np.meshgrid(np.arange(x0, x0 + Nx, dtype=np.int64),
                         np.arange(y0, y0 + Ny, dtype=np.int64))
    good_pix = ~np.isnan(img)
    
    return coo_matrix((img[good_pix], (iy[good_pix], ix[good_pix])),
                      shape=(Ny_tot, Nx_tot))
