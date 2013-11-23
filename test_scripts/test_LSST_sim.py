"""
Testing WCS projections on LSST simulation files
"""
import os, sys
sys.path.append(os.path.abspath('..'))


import numpy as np
import matplotlib.pyplot as plt

from spheredb.get_data import\
    get_stripe82_file, all_lsst_exposures, get_LSST_file
from spheredb.conversions import FITS_to_HPX, HPX_grid_step
from spheredb.util import regrid

import os
import pyfits

import re
import datetime

# Note: USE INSERT NOT MERGE!!!!


if 1:
    from scidbpy import interface
    sdb = interface.SciDBShimInterface('http://vega.cs.washington.edu:8080')
    Nside = 2 ** 16 #19
    hdulist = get_LSST_file()
    output = FITS_to_HPX(hdulist[1].header, hdulist[1].data, Nside,
                         return_sparse=True)

    print output.shape

    RA_range = (output.row.min(), output.row.max())
    DEC_range = (output.col.min(), output.col.max())

    dRA = RA_range[1] - RA_range[0]
    dDEC = DEC_range[1] - DEC_range[0]

    RA_range = (RA_range[0] - 1 * dRA, RA_range[1] + 1 * dRA)
    DEC_range = (DEC_range[0] - 1 * dDEC, DEC_range[1] + 1 * dDEC)

    arr = sdb.from_sparse(output)

    subarr = arr[RA_range[0]:RA_range[1],
                 DEC_range[0]:DEC_range[1]]

    plt.imshow(np.log(subarr.toarray()), cmap=plt.cm.binary)
    plt.show()

elif 1:
    times = [hdulist[1].header['TAI'] for hdulist in all_lsst_exposures()]
    times = np.asarray(times)
    times.sort()
    print times.min()
    print times.max()
    plt.plot(24 * (times - 50095), '.k')
    plt.show()
