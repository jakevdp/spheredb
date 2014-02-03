import numpy as np

from .lsst_warp import LSSTWarper
from scidbpy import interface

SHIM_DEFAULT = 'http://localhost:8080'


def find_index_bounds(arr, sdb, dims=None):
    """Find the bounds of the nonempty data in an array

    Parameters
    ----------
    arr : SciDBArray
        the input scidb array
    sdb : SciDBInterface
        interface to scidb server
    dims : int, tuple, or None
    """
    shape = arr.shape
    if dims is None:
        dims = range(len(shape))
    else:
        try:
            dims = tuple(dims)
        except TypeError:
            dims = (dims,)

    if(min(dims) < 0 or max(dims) >= len(shape)):
        raise ValueError("dims out of range")

    # todo: we need some unique dimension names here
    #       there should be a scidbpy utility for this
    dim_names = ["_tmp{0:01d}".format(d) for d in dims]
    query_string = "store(aggregate(apply({A},"
    query_string += ', '.join("{0}, {{A.d{1}}}".format(d, i)
                             for (d, i) in zip(dim_names, dims))
    query_string += "),"
    query_string += ", ".join("min({0}), max({0})".format(d)
                              for d in dim_names)
    query_string += "), {output})"

    output = sdb.new_array()
    sdb.query(query_string, A=arr, output=output)

    output = output.toarray()
    return np.asarray([output[name][0] for name in output.dtype.names])


class HPXPixels3D(object):
    """
    Class to store and interact with 3D Healpix-projected data

    The three dimensions include two angular dimensions and one time dimension.
    """
    def __init__(self, name=None, input_files=None,
                 cdelt=3, cunit='arcsec', kernel='lanczos2',
                 force_reload=False, interface=None):
        self.name = name
        self.force_reload = force_reload
        self.interface = interface

        if self.interface is None:
            self.interface = self.open_scidb_connection()

        self.warper = LSSTWarper(cdelt=cdelt,
                                 cunit=cunit,
                                 kernel=kernel,
                                 interface=self.interface)

        if (name is not None):
            arr_exists = (name in self.interface.list_arrays())
                
            if force_reload or not arr_exists:
                print "loading into array: {0}".format(self.name)

                # clear the array if needed
                if arr_exists:
                    self.interface.query("remove({0})", name)

                self.arr = None
                self._load_files(input_files)
            else:
                print "using existing array: {0}".format(self.name)
                self.arr = self.interface.wrap_array(self.name)
        else:
            self.arr = None
            self._load_files(input_files)

    @staticmethod
    def open_scidb_connection(address=SHIM_DEFAULT):
        return interface.SciDBShimInterface(address)

    def _load_files(self, files):
        for i, fitsfile in enumerate(files):
            print "- ({0}/{1}) loading {2}".format(i + 1,
                                                   len(files),
                                                   fitsfile)
                                                   
            if self.arr is None:
                self.arr = self.warper.scidb3d_from_fits(fitsfile)
                if self.name is not None:
                    self.arr.rename(self.name, persistent=True)
            else:
                self.interface.query("insert({0}, {1})",
                                     self.warper.scidb3d_from_fits(fitsfile),
                                     self.arr)

    def time_slice(self, time1, time2=None):
        if time2 is None:
            return HPXPixels2D(self, self.arr[:, :, time1])
        else:
            raise NotImplementedError()

    def coadd(self):
        return HPXPixels2D(self, self.arr.sum(2))

    def unique_times(self):
        return self.arr.max((0, 1)).tosparse()['time']

    def index_bounds(self):
        bounds = find_index_bounds(self.arr, self.interface)
        return bounds[:2], bounds[2:4], bounds[4:6]

class HPXPixels2D(object):
    """Container for 2D LSST Pixels stored in SciDB"""
    def __init__(self, pix3d, arr):
        self.pix3d = pix3d
        self.arr = arr

    @property
    def interface(self):
        return self.pix3d.interface

    def regrid(self, *args):
        """Parameters
        
        N: number of pixels to add in each dimesion
        aggregrate: scidb aggregate to use
        """
        return HPXPixels2D(self.pix3d, self.arr.regrid(*args))

    def subarray(self, xlim, ylim):
        return HPXPixels2D(self.pix3d, self.arr[xlim[0]:xlim[1],
                                                ylim[0]:ylim[1]])

    def index_bounds(self):
        bounds = find_index_bounds(self.arr, self.interface)
        return bounds[:2], bounds[2:4]
