import os, sys
sys.path.append(os.path.abspath('..'))

from spheredb.sdb_from_warped import sdb_from_warped

arr = sdb_from_warped('out.fits')
print arr.shape
print arr.nnz
