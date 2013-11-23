"""
LSST Warping

This requires the LSST stack to be installed and setup: see
https://dev.lsstcorp.org/trac/wiki/Installing

After install, run the following (adapting for your install path):

source ~/LSST_STACK/loadLSST.sh
setup python
setup afw
"""
import os, sys
sys.path.append(os.path.abspath('..'))

import numpy as np
#from scipy.sparse import coo_matrix
from lsst import afw
import lsst.afw.image as afwImage
import lsst.afw.math as afwMath
import lsst.daf.base as dafBase

from spheredb.lsst_warp import LSSTWarper

W = LSSTWarper()
filename = "/Users/jakevdp/research/LSST_IMGS/v865833781-fr/R21/S12.fits"

W.warp_and_save(infile=filename, outfile='out.fits')
exit()


def make_HPX_wcs(cdelt, cunit='deg'):
    ps = dafBase.PropertySet()
    ps.add('NAXIS',2)
    ps.add('CTYPE1', 'RA---HPX')
    ps.add('CTYPE2', 'DEC--HPX')
    ps.add('CUNIT1', cunit)
    ps.add('CUNIT2', cunit)
    ps.add('CDELT1', cdelt)
    ps.add('CDELT2', cdelt)
    ps.add('CRVAL1', 0)
    ps.add('CRVAL2', 0)
    ps.add('CRPIX1', 0)
    ps.add('CRPIX2', 0)
    return afwImage.makeWcs(ps)


def read_fits_and_warp(fitsfile, cdelt=1, cunit='arcsec',
                       kernel='lanczos2'):
    exp = afwImage.ExposureF(fitsfile)
    wcs_in = exp.getWcs()
    wcs_out = make_HPX_wcs(cdelt, cunit)
    warper = afwMath.Warper(kernel)
    warpedExposure = warper.warpExposure(destWcs=wcs_out,
                                         srcExposure=exp)
    return warpedExposure


def warped_to_sparse(fitsfile, cdelt=0.001, kernel='lanczos2'):
    # hard-code cunit for our indices
    cunit = 'deg'

    warped = read_fits_and_warp(fitsfile, cdelt=cdelt,
                                cunit=cunit, kernel=kernel)
    img = warped.getMaskedImage()
    x0, y0 = img.getXY0()
    y0 += int(np.round(90 / cdelt))

    Nx = img.getWidth()
    Nx_tot = int(np.round(180. / cdelt))

    Ny = img.getHeight()
    Ny_tot = int(np.round(90. / cdelt))

    img, mask, err = img.getArrays()


    ix = np.arange(x0, x0 + Nx, dtype=np.int64)
    iy = np.arange(y0, y0 + Ny, dtype=np.int64)
    ix, iy = np.meshgrid(ix, iy)

    ix, iy, img = map(np.ravel, (ix, iy, img))
    good_pixels = ~np.isnan(img)
    ix = ix[good_pixels]
    iy = iy[good_pixels]
    img = img[good_pixels]
    
    return coo_matrix((im, (iy, ix)),
                      shape=(Ny, Nx))


if __name__ == '__main__':
    filename = "/Users/jakevdp/research/LSST_IMGS/v865833781-fr/R21/S12.fits"
    sp = warped_to_sparse(filename)

    print sp.shape
    print sp.nnz

    #warped = read_fits_and_warp(filename)
    #warped.writeFits('out.fits')
    #import IPython; IPython.embed()
