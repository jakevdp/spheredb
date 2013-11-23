"""Tools to fetch stripe 82 data"""

import os
import re
import urllib2
from bz2 import BZ2File
from StringIO import StringIO

import pyfits
from astroML.datasets.tools import download_with_progress_bar


URL = "http://data.sdss3.org/dr10/env/BOSS_PHOTOOBJ/frames/{rerun}/{run}/{camcol}/frame-{filter}-{run:06d}-{camcol}-{framenum:04d}.fits.bz2"

LSST_DIR = '../../LSST_IMGS'


def get_stripe82_url(rerun, run, camcol=1, filter='u', framenum=1):
    """Return a Stripe 82 FITS URL"""
    return URL.format(run=run, rerun=rerun, camcol=camcol,
                      filter=filter, framenum=framenum)


def get_stripe82_file(rerun, run, camcol=1, filter='u', framenum=1):
    """Get a Stripe 82 FITS file"""
    data_url = get_stripe82_url(rerun, run, camcol, filter, framenum)
    local_file = os.path.split(data_url)[-1]

    if os.path.exists(local_file):
        print "using local image:", local_file
    else:
        buf = download_with_progress_bar(data_url, return_buffer=True)
        open(local_file, 'wb').write(buf.read())

    bz2f = BZ2File(local_file)
    hdulist = pyfits.open(StringIO(bz2f.read()))

    return hdulist


def all_lsst_exposures(lsst_dir=LSST_DIR):
    R = re.compile("^S[0-2][0-2]\.fits$")
    exposures = []
    for dirpath, dirnames, filenames in os.walk(lsst_dir):
        for f in filter(R.match, filenames):
            hdulist = pyfits.open(os.path.join(dirpath, f))
            yield hdulist



def get_LSST_file(lsst_dir=LSST_DIR,
                  loc='v865833781-fr/R21/S12.fits'):
    R = re.compile("^S[0-2][0-2]\.fits$")

    path = os.path.join(os.path.abspath(lsst_dir), loc)
    if not os.path.exists(path):
        raise ValueError("Cannot find file {0}".format(path))
    hdulist = pyfits.open(path)
    return hdulist


def _plot_test():
    import numpy as np
    import matplotlib.pyplot as plt

    framenum=18

    for rerun in [2700, 2703]:
        hdulist = get_file(301, rerun, 1, 'g', framenum)
        im = hdulist[0].data
        plt.figure()
        plt.imshow(np.log10(im), cmap=plt.cm.binary)
    
    plt.show()
