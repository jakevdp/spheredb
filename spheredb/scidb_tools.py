from scidbpy import interface
from .conversions import FITS_to_HPX, HPX_grid_size

SHIM_DEFAULT = 'http://localhost:8080'
NSIDE_DEFAULT = 2 ** 16

# TODO: these should be derived from SciDBArray (carefully...)


class LSSTPixels3D(object):
    """Container for 3D LSST Pixels stored in SciDB"""
    def __init__(self, name,
                 files=None, hdu_index=1, Nside=NSIDE_DEFAULT,
                 sdb=None, quiet=False, force_reload=False):
        if sdb is None:
            self.sdb = self.open_scidb_connection()
        else:
            self.sdb = sdb

        self.quiet = quiet
        self.Nside = Nside
        self.arr = None

        if ((not force_reload) and (name is not None)
                  and (name in self.sdb.list_arrays())):
            print "Using existing array {0}".format(name)
            self.arr = self.sdb.wrap_array(name)

        elif files is None:
            raise ValueError("must specify either existing name, "
                             "or pass file iterator")

        else:
            print "Loading array from files"
            self._load_files(files, hdu_index)

    def _load_files(self, files, hdu_index):
        Nx_hpx, Ny_hpx = HPX_grid_size(self.Nside)
        Ntimes = int(100000 * 24 * 60 * 60)  # essentially unlimited
        shape = (Nx_hpx, Ny_hpx, Ntimes)
        dim_names = ('x', 'y', 'time')
        dtype = '<val:double>'

        for hdulist in files:
            if not quiet:
                print "uploading file at", hdulist[1].header['TAI']
            hdu = hdulist[hdu_index]
            output = FITS_to_HPX(hdu.header, hdu.data,
                                 Nside, return_sparse=False)
                             
            new_arr = self.sdb.from_array(output)
            redimensioned_arr = self.sdb.new_array(shape, dtype,
                                                   dim_names=dim_names)
            self.sdb.query("redimension_store({0}, {1})",
                           new_arr, redimensioned_arr)

            # First time around
            if self.arr is None:
                self.arr = redimensioned_arr
                self.arr.rename(name, persistent=True)
            else:
                self.sdb.query("insert({0}, {1})",
                               redimensioned_arr, self.arr)

    def unique_times(self):
        return self.arr.max((0, 1)).tosparse()['time']

    def time_slice(self, time1, time2=None):
        if time2 is None:
            return LSSTPixels2D(self, self.arr[:, :, time1])
        else:
            raise NotImplementedError()

    def coadd(self):
        return LSSTPixels2D(self, self.arr.sum(2))

    @staticmethod
    def open_scidb_connection(address=SHIM_DEFAULT):
        return interface.SciDBShimInterface(address)


class LSSTPixels2D(object):
    """Container for 2D LSST Pixels stored in SciDB"""
    def __init__(self, pix3d, arr):
        self.pix3d = pix3d
        self.sdb = pix3d.sdb
        self.arr = arr

    def regrid(self, *args):
        return LSSTPixels2D(self.pix3d, self.arr.regrid(*args))

    def subarray(self, xlim, ylim):
        return LSSTPixels2D(self.pix3d, self.arr[xlim[0]:xlim[1],
                                                 ylim[0]:ylim[1]])
