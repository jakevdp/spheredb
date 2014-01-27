"""
LSST Warping

This requires the LSST stack to be installed and setup: see
https://dev.lsstcorp.org/trac/wiki/Installing

After install, run the following (adapting for your install path):

[~]$ source ~/LSST_STACK/loadLSST.sh
[~]$ setup python
[~]$ setup afw
"""
import lsst.afw.image as afwImage
import lsst.afw.math as afwMath
import lsst.daf.base as dafBase


class LSSTWarper(object):
    def __init__(self, cunit='arcsec', cdelt=1, kernel='lanczos2'):
        self.kernel = kernel
        cunit = cunit.lower().strip()
        if cunit == 'deg':
            self.cdelt = cdelt
        elif cunit == 'arcmin':
            self.cdelt = cdelt * 1. / 60.
        elif cunit == 'arcsec':
            self.cdelt = cdelt * 1. / 3600.
        else:
            raise ValueError("Unrecognized cunit: {0}".format(cunit))

    @staticmethod
    def make_wcs(cdelt, cunit='deg'):
        """Construct a HEALPix WCS header"""
        ps = dafBase.PropertySet()
        ps.add('NAXIS', 2)
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

    def warped_from_fits(self, fitsfile):
        """Return a warped exposure computed from an LSST exposure"""
        exp = afwImage.ExposureF(fitsfile)
        wcs_in = exp.getWcs()
        wcs_out = self.make_wcs(self.cdelt, 'deg')
        warper = afwMath.Warper(self.kernel)
        warpedExposure = warper.warpExposure(destWcs=wcs_out,
                                             srcExposure=exp)
        return warpedExposure

    def warp_and_save(self, infile, outfile):
        warpedExposure = self.warped_from_fits(infile)
        warpedExposure.writeFits(outfile)

    def sparse_from_fits(self, fitsfile):
        """Return a sparse HPX array from an LSST exposure"""
        import numpy as np
        from scipy.sparse import coo_matrix
        warped = self.warped_from_fits(fitsfile)

        img = warped.getMaskedImage()
        x0, y0 = img.getXY0()
        y0 += int(np.round(90 / self.cdelt))

        Nx = img.getWidth()
        Nx_tot = int(np.round(180. / self.cdelt))

        Ny = img.getHeight()
        Ny_tot = int(np.round(90. / self.cdelt))

        img, mask, err = img.getArrays()

        ix = np.arange(x0, x0 + Nx, dtype=np.int64)
        iy = np.arange(y0, y0 + Ny, dtype=np.int64)
        ix, iy = np.meshgrid(ix, iy)

        ix, iy, img = map(np.ravel, (ix, iy, img))
        good_pixels = ~np.isnan(img)
        ix = ix[good_pixels]
        iy = iy[good_pixels]
        img = img[good_pixels]

        return coo_matrix((img, (iy, ix)),
                          shape=(Ny_tot, Nx_tot))
