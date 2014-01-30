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
HPX_data = HPXPixels3D(input_files=filenames[:2])
print HPX_data.arr.shape

time = HPX_data.unique_times()[0]
sliced = HPX_data.time_slice(time)
print sliced.arr.shape

coadded = HPX_data.coadd()
print coadded.arr.shape
