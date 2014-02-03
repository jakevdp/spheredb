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

import matplotlib.pyplot as plt
import numpy as np

sys.path.append(os.path.abspath('..'))
from spheredb.scidb_tools import HPXPixels3D, find_index_bounds

filenames = glob.glob("/home/jakevdp/research/LSST_IMGS/*/R*/S*.fits")
print "total number of files:", len(filenames)

HPX_data = HPXPixels3D(input_files=filenames[:20],
                       name='LSSTdata', force_reload=False)
times = HPX_data.unique_times()

xlim, ylim, tlim = HPX_data.index_bounds()


for time in times[:2]:
    tslice = HPX_data.time_slice(time)
    tslice_arr = tslice.arr[xlim[0]:xlim[1],
                            ylim[0]:ylim[1]].toarray()
    fig, ax = plt.subplots()
    im = ax.imshow(np.log(tslice_arr), cmap=plt.cm.binary)
    ax.set_xlim(400, 440)
    ax.set_ylim(860, 820)
    fig.colorbar(im, ax=ax)
    ax.set_title("time = {0}".format(time))

coadd = HPX_data.coadd().arr[xlim[0]:xlim[1],
                             ylim[0]:ylim[1]].toarray()
fig, ax = plt.subplots()
im = ax.imshow(np.log(coadd), cmap=plt.cm.binary)
ax.set_xlim(400, 440)
ax.set_ylim(860, 820)
fig.colorbar(im, ax=ax)
ax.set_title("coadd")


plt.show()
