from .lsst_warp import LSSTWarper
from scidbpy import interface

SHIM_DEFAULT = 'http://localhost:8080'


class HPXPixels3D(object):
    """
    Class to store and interact with 3D Healpix-projected data

    The three dimensions include two angular dimensions and one time dimension.
    """
    def __init__(self, name=None, input_files=None,
                 cdelt=3, cunit='arcsec', kernel='lanczos2',
                 force_reload=False, sdb=None):
        self.name = name
        self.warper = LSSTWarper(cdelt=cdelt,
                                 cunit=cunit,
                                 kernel=kernel)
        self.force_reload = force_reload
        self.sdb = sdb

        if self.sdb is None:
            self.sdb = self.open_scidb_connection()

        # create an empty array here?
        self.arr = None

        self._load_files(input_files)
            
    @staticmethod
    def open_scidb_connection(address=SHIM_DEFAULT):
        return interface.SciDBShimInterface(address)

    def _load_files(self, files):
        for fitsfile in files:
            arr = self.warper.scidb3d_from_fits(fitsfile, self.sdb)

            if self.arr is None:
                self.arr = arr
            else:
                self.sdb.query("insert({0}, {1})",
                               arr, self.arr)

    def unique_times(self):
        return self.arr.max((0, 1)).tosparse()['time']
