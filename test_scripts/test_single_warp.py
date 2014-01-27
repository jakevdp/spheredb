"""
Single LSST Warp
----------------
This tests the warping of a single LSST exposure into a sparse matrix
representation of a HealPix grid.
"""
import os, sys
sys.path.append(os.path.abspath('..'))
from spheredb.lsst_warp import LSSTWarper

filename = "~/research/LSST_IMGS/v865833781-fr/R21/S12.fits"

print "loading frame from {0}".format(filename)
W = LSSTWarper(cdelt=3, cunit='arcsec')
sp = W.sparse_from_fits(filename)

print "Shape of HPX sparse array:", sp.shape
print "number of nonzero entries:", sp.nnz
