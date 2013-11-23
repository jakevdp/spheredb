"""
LSST Warping

This requires the LSST stack to be installed and setup: see
https://dev.lsstcorp.org/trac/wiki/Installing

After install, run the following (adapting for your install path):

source ~/LSST_STACK/loadLSST.sh
setup python
setup afw
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
