"""
Plot LSST Frame
---------------
Load an LSST frame into scidb and plot the results
"""
import os, sys
sys.path.append(os.path.abspath('..'))

# 1. Use spheredb to load the LSST exposure and warp onto HPX
from spheredb.lsst_warp import LSSTWarper
filename = "~/research/LSST_IMGS/v865833781-fr/R21/S12.fits"
print "loading frame from {0}".format(filename)
W = LSSTWarper(cdelt=3, cunit='arcsec')
sp = W.sparse_from_fits(filename)

# 2. Load into SciDB to regrid
from scidbpy import interface
sdb = interface.SciDBShimInterface('http://localhost:8080')
M = sdb.from_sparse(sp)
R = M.regrid(1000, aggregate="sum")

# 3. Use matplotlib to plot the down-sampled version
import numpy as np
import matplotlib.pyplot as plt

plt.imshow(R.toarray(), interpolation='nearest', cmap=plt.cm.binary)
plt.show()

