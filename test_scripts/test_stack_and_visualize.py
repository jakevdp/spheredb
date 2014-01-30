"""
Stacking and Visualizing
------------------------
This script does the following:

1. Input LSST images, warp to sparse matrix, store as scidb arrays.
This tests the warping of a single LSST exposure into a sparse matrix
representation of a HEALPix grid.
"""
import os
import sys
import glob
sys.path.append(os.path.abspath('..'))
from spheredb.scidb_tools import HPXPixels3D

filenames = glob.glob("/home/jakevdp/research/LSST_IMGS/*/R*/S*.fits")
print "total number of files:", len(filenames)

HPX_data = HPXPixels3D(input_files=filenames[:20],
                       name='LSSTdata', force_reload=False)
times = HPX_data.unique_times()

for time in times:
    tslice = HPX_data.time_slice(time)
    print time, tslice.arr.shape, tslice.arr.nonempty()
