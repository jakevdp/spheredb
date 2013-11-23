__all__ = ['HPX_grid_step', 'HPX_grid_size', 'FITS_to_HPX']

import numpy as np
from scipy import sparse

# Kapteyn software contains tie-ins to WCS standard.
try:
    from kapteyn import wcs
except ImportError:
    print ("kapteyn package required: download at\n"
           "http://www.astro.rug.nl/software/kapteyn/")
    raise

from .grid_interpolation import GridInterpolation
from .util import regrid


def HPX_grid_size(Nside):
    """Return the size of the pixel grid (Nx, Ny) for a given Nside"""
    Nx = 8 * Nside
    Ny = 4 * Nside + 1
    return Nx, Ny


def HPX_grid_step(Nside):
    """Return the size of the step between pixels in degrees"""
    return 45. / Nside


def FITS_to_HPX(header, data, Nside, return_sparse=False):
    """Convert data from FITS format to sparse HPX grid

    Parameters
    ----------
    header : dict or PyFITS header
        WCS header describing the coordinates of the input array
    data : array_like
        Input data array
    Nside : int
        HEALPix gridding parameter

    Returns
    -------
    hpx_data : csr matrix
        The HPX-projected data
    """
    # Here's what we do for this function: we're working in "IMG coords"
    # (i.e. the projection of the input data) and "HPX coords" (i.e. the
    # projection of the output data).  In between, we use "WCS coords".
    #
    # These are the steps involved:
    #  1. Create an array of image edge-pixels in IMG coords, and project
    #     these to HPX coords.
    #  2. From these bounds, create a regular grid of HPX coords that covers
    #     the image.  Project this grid to IMG coords.
    #  3. In IMG coords, interpolate the image data to the healpix grid.
    #  4. Use this data to construct a sparse array in HPX coords.

    if header['NAXIS'] != 2:
        raise ValueError("input data & header must be two dimensional")

    if data.shape != (header['NAXIS2'], header['NAXIS1']):
        raise ValueError("data shape must match header metadata")

    # Create wcs projection instance from the header
    proj_img = wcs.Projection(header)
    
    # Create wcs projection for healpix grid
    # Note that the "pixel" coordinates here are measured in degrees...
    # 0 to 360 in x/RA and -90 to 90 in y/DEC
    proj_hpx = wcs.Projection({'NAXIS': 2,
                               'CTYPE1': 'RA---HPX',
                               'CTYPE2': 'DEC--HPX'})

    # Define the dimension of the HEALPIX SciDB grid
    Nx_hpx, Ny_hpx = HPX_grid_size(Nside)
    dx_hpx = dy_hpx = HPX_grid_step(Nside)
    
    #x_hpx = np.linspace(0, 360, Nx_hpx, endpoint=False)
    #y_hpx = np.linspace(-90, 90, Ny_hpx)

    # Find the coordinates of the pixels at the edge of the image
    # Projecting these onto the healpix grid will give the bounds we need.
    img_bounds_x = np.arange(header['NAXIS2'])
    zeros_x = np.zeros_like(img_bounds_x)
    img_bounds_y = np.arange(header['NAXIS1'])
    zeros_y = np.zeros_like(img_bounds_y)

    img_bounds_pix = np.concatenate(
        [img_bounds_x, img_bounds_x, zeros_y, zeros_y + img_bounds_x[-1],
         zeros_x, zeros_x + img_bounds_y[-1], img_bounds_y, img_bounds_y]
    ).reshape((2, -1)).T

    x_bound_hpx, y_bound_hpx =\
                    proj_hpx.topixel(proj_img.toworld(img_bounds_pix)).T    

    # here we take the pixels at the edge of the boundaries of the image,
    # transform them to HPX coordinates, and find the required extent
    # of the HPX pixel grid.
    #    [TODO: check for crossing the pole]
    # first we need to calculate pixel number
    i_bound_hpx = x_bound_hpx / dx_hpx
    j_bound_hpx = (y_bound_hpx + 90.) / dy_hpx
    
    i_hpx = np.arange(int(np.floor(i_bound_hpx.min())),
                      int(np.ceil(i_bound_hpx.max()) + 1))
    j_hpx = np.arange(int(np.floor(j_bound_hpx.min())),
                      int(np.ceil(j_bound_hpx.max()) + 1))

    x_hpx = i_hpx * dx_hpx
    y_hpx = j_hpx * dy_hpx - 90.
    
    # Create the grid of HPX pixels
    pixel_ind_hpx = np.vstack(map(np.ravel, np.meshgrid(i_hpx, j_hpx))).T
    pixel_locs_hpx = np.vstack(map(np.ravel, np.meshgrid(x_hpx, y_hpx))).T
    pixel_locs_img = proj_img.topixel(proj_hpx.toworld(pixel_locs_hpx))

    ## DEBUG: Plot the borders & grid in the HPX projection
    #import matplotlib.pyplot as plt
    #plt.plot(i_bound_hpx, j_bound_hpx, '.k')
    #plt.plot(pixel_ind_hpx[:, 0], pixel_ind_hpx[:, 1], '.r')
    #plt.show()
    #exit()

    ## DEBUG: Plot the HPX grid in the IMG projection
    #import matplotlib.pyplot as plt
    #plt.plot(img_bounds_pix[:, 0], img_bounds_pix[:, 1], '.k')
    #plt.plot(pixel_locs_img[:, 0], pixel_locs_img[:, 1], '.r')
    #plt.show()
    #exit()

    # Interpolate from data to pixel locations
    I = GridInterpolation(data, [0, 0], [1, 1])
    HPX_vals = I(pixel_locs_img)#.reshape(len(y_hpx), len(x_hpx))

    # # DEBUG: Plot regridded input data next to the interpolated HPX data
    # import matplotlib.pyplot as plt
    # plt.figure(figsize=(8, 8))
    # plt.subplot(211, aspect='equal')
    # plt.contourf(x_hpx, y_hpx, HPX_vals)
    # plt.subplot(212, aspect='equal')
    # plt.contourf(regrid(data, 5))
    # plt.show()
    # exit()
    
    good_vals = ~np.isnan(HPX_vals)
    x, y = pixel_ind_hpx[good_vals].T
    HPX_vals = HPX_vals[good_vals]

    if return_sparse:
        return sparse.coo_matrix((HPX_vals, (x, y)),
                                 shape=(Nx_hpx, Ny_hpx))
    else:
        output = np.zeros(len(HPX_vals),
                          dtype=[('time', np.int64),
                                 ('x', np.int64),
                                 ('y', np.int64),
                                 ('val', np.float64)])
        # use MJD in seconds
        output['time'] = int(header['TAI'] * 24 * 60 * 60)
        output['x'] = x
        output['y'] = y
        output['val'] = HPX_vals
        return output
