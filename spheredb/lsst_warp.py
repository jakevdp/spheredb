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

import numpy as np


class LSSTWarper(object):
    """Tools to warp input fits data to a HEALPix grid."""
    def __init__(self, cunit='arcsec', cdelt=1, kernel='lanczos2',
                 interface=None):
        self.kernel = kernel
        self.cdelt = cdelt
        self.cunit = cunit.lower().strip()
        self.interface = interface
        if self.cunit not in ['deg', 'arcmin', 'arcsec']:
            raise ValueError("cunit='{0}' not recognized".format(self.cunit))

    @classmethod
    def compute_cdelt_deg(cls, cdelt, cunit):
        cunit = cunit.lower().strip()
        if cunit == 'deg':
            return cdelt
        elif cunit == 'arcmin':
            return cdelt * 1. / 60.
        elif cunit == 'arcsec':
            return cdelt * 1. / 3600.
        else:
            raise ValueError("Unrecognized cunit: {0}".format(cunit))

    @classmethod
    def grid_size(cls, cdelt, cunit):
        """Return Nx, Ny for the given cdelt and cunit"""
        cdelt = cls.compute_cdelt_deg(cdelt, cunit)
        Nx = int(np.round(180. / cdelt))
        Ny = int(np.round(90. / cdelt))
        return (Nx, Ny)

    @property
    def cdelt_deg(self):
        return self.compute_cdelt_deg(self.cdelt, self.cunit)

    @property
    def Nx(self):
        return int(np.round(180. / self.cdelt_deg))

    @property
    def Ny(self):
        return int(np.round(90. / self.cdelt_deg))

    @property
    def Nt(self):
        return int(100000 * 24 * 60 * 60)

    def make_wcs(self):
        """Construct a HEALPix WCS header"""
        ps = dafBase.PropertySet()
        ps.add('NAXIS', 2)
        ps.add('CTYPE1', 'RA---HPX')
        ps.add('CTYPE2', 'DEC--HPX')
        ps.add('CUNIT1', 'deg')
        ps.add('CUNIT2', 'deg')
        ps.add('CDELT1', self.cdelt_deg)
        ps.add('CDELT2', self.cdelt_deg)
        ps.add('CRVAL1', 0)
        ps.add('CRVAL2', 0)
        ps.add('CRPIX1', 0)
        ps.add('CRPIX2', 0)
        return afwImage.makeWcs(ps)

    def get_exposure_date(self, fitsfile):
        metadata = afwImage.ExposureF(fitsfile).getMetadata()
        return metadata.get('MJD-OBS')

    def warped_from_fits(self, fitsfile):
        """Return a warped exposure computed from an LSST exposure"""
        exp = afwImage.ExposureF(fitsfile)
        wcs_in = exp.getWcs()
        wcs_out = self.make_wcs()
        warper = afwMath.Warper(self.kernel)
        warpedExposure = warper.warpExposure(destWcs=wcs_out,
                                             srcExposure=exp)
        return warpedExposure

    def warp_and_save(self, infile, outfile):
        warpedExposure = self.warped_from_fits(infile)
        warpedExposure.writeFits(outfile)

    def sparse_from_fits(self, fitsfile):
        """Return a sparse HPX array from an LSST exposure"""
        from scipy import sparse

        warped = self.warped_from_fits(fitsfile)

        img = warped.getMaskedImage()
        x0, y0 = img.getXY0()
        y0 += self.Ny

        Nx_img = img.getWidth()
        Ny_img = img.getHeight()

        img, mask, err = img.getArrays()

        ix = np.arange(x0, x0 + Nx_img, dtype=np.int64)
        iy = np.arange(y0, y0 + Ny_img, dtype=np.int64)
        ix, iy = np.meshgrid(ix, iy)

        ix, iy, img = map(np.ravel, (ix, iy, img))
        good_pixels = ~np.isnan(img)
        ix = ix[good_pixels]
        iy = iy[good_pixels]
        img = img[good_pixels]

        return sparse.coo_matrix((img, (iy, ix)),
                                 shape=(self.Ny, self.Nx))

    def scidb2d_from_fits(self, filename):
        """Return a SciDB array from a fits file"""
        if self.interface is None:
            raise ValueError("scidb interface must be defined")

        sp = self.sparse_from_fits(filename)
        return self.interface.from_sparse(sp)

    def scidb3d_from_fits(self, fitsfile):
        if self.interface is None:
            raise ValueError("scidb interface must be defined")
            
        time = self.get_exposure_date(fitsfile)
        warped = self.sparse_from_fits(fitsfile)
            
        warped_data = np.zeros(warped.nnz, dtype=[('time', np.int64),
                                                  ('x', np.int64),
                                                  ('y', np.int64),
                                                  ('val', np.float64)])

        warped_data['time'] = int(time * 24 * 60 * 60)
        warped_data['x'] = warped.row
        warped_data['y'] = warped.col
        warped_data['val'] = warped.data

        warped_arr = self.interface.from_array(warped_data)
        redimensioned = self.interface.new_array(shape=(self.Nx, self.Ny,
                                                        self.Nt),
                                                 dtype='<val:double>',
                                                 dim_names=('x', 'y', 'time'))
        self.interface.query('redimension_store({0}, {1})',
                             self.interface.from_array(warped_data),
                             redimensioned)
        return redimensioned
        
