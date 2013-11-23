"""
Warp all the LSST example exposures, saving them as warped FITS files.

This requires the LSST stack to be installed and setup: see
https://dev.lsstcorp.org/trac/wiki/Installing

After install, run the following (adapting for your install path):

source ~/LSST_STACK/loadLSST.sh
setup python
setup afw
"""
import os, sys
sys.path.append(os.path.abspath('..'))
import re

from spheredb.lsst_warp import LSSTWarper


#from spheredb.get_data import all_lsst_files
LSST_DIR = os.path.expanduser("~/research/LSST_IMGS")
def all_lsst_files(lsst_dir=LSST_DIR):
    R = re.compile("^S[0-2][0-2]\.fits$")
    exposures = []
    for dirpath, dirnames, filenames in os.walk(lsst_dir):
        for f in filter(R.match, filenames):
            yield os.path.join(dirpath, f)


def main():
    if not os.path.exists('warped'):
        os.makedirs('warped')

    W = LSSTWarper()
    for i, f in enumerate(all_lsst_files()):
        print " - warping {0}".format(f)
        W.warp_and_save(infile=f,
                        outfile='warped/out{0:03d}.fits'.format(i + 1))


if __name__ == '__main__':
    main()
