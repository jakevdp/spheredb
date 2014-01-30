import numpy as np
from numpy.testing import assert_equal, assert_raises
from scipy import sparse

from scidbpy import interface
from spheredb.scidb_tools import SHIM_DEFAULT, find_index_bounds

sdb = interface.SciDBShimInterface(SHIM_DEFAULT)

def test_index_bounds():
    row = [4, 4, 5, 5]
    col = [3, 6, 2, 4]
    data = [1, 1, 1, 1]
    M = sparse.coo_matrix((data, (row, col)), shape=(10, 10))
    Msdb = sdb.from_sparse(M)

    assert_equal(find_index_bounds(Msdb, sdb), [4, 5, 2, 6])
    assert_equal(find_index_bounds(Msdb, sdb, 1), [2, 6])
    assert_equal(find_index_bounds(Msdb, sdb, (0, 1)), [4, 5, 2, 6])
    assert_raises(ValueError, find_index_bounds, Msdb, sdb, 3)
