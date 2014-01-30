"""
Plot LSST Frame
---------------
Load an LSST frame into scidb and plot the results
"""
import os, sys
sys.path.append(os.path.abspath('..'))

# 1. Set up LSST Warper
from spheredb.lsst_warp import LSSTWarper
filename = "~/research/LSST_IMGS/v865833781-fr/R21/S12.fits"
print "loading frame from {0}".format(filename)
W = LSSTWarper(cdelt=3, cunit='arcsec')

# 2. Push Into SciDB
from scidbpy import interface
sdb = interface.SciDBShimInterface('http://localhost:8080')
M = W.scidb_from_fits(filename, sdb)
R = M.regrid(1000, aggregate="sum")

# 3. Use matplotlib to plot the down-sampled version
import numpy as np
import matplotlib.pyplot as plt

plt.imshow(R.toarray(), interpolation='nearest', cmap=plt.cm.binary)
plt.show()

