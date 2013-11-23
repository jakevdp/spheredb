import numpy as np
import itertools


def coo_to_recarray(M):
    dtype = [('data', np.float), ('i1', np.int), ('i2', np.int)]
    A = np.empty(len(M.data), dtype=dtype)
    A['data'] = M.data
    A['i1'] = M.row
    A['i2'] = M.col

    return A


def regrid(X, D, agg=np.sum):
    """Regrid an N-dimensional matrix using numpy aggregates

    Parameters
    ----------
    X : array_like
        N-dimensional data to regrid
    D : int or array-like
        The step size of the regridding.  D must be either an integer
        or an array of integers of length X.ndim
    agg : function (optional)
        numpy aggregate operation (default=np.sum). Must be a function which
        takes an array as the first argument, and dimensions to aggregate over
        as the second argument.

    Returns
    -------
    X_agg : numpy array
        for X.shape = (N1, N2, N3...) and D = (D1, D2, D3...) then the shape
        of X_agg is (N1 // D1, N2 // D2, N3 // D3...)

    Note
    ----
    If Ni is not a multiple of Di, then it will be truncated from the right.
    """
    X = np.asarray(X)
    D = np.zeros(X.ndim, dtype=int) + np.asarray(D).astype(int)
    N = np.asarray(X.shape)
    final_shape = N // D
    input_shape = D * (N // D)

    # Truncate the edges of X
    X = X[[slice(None, s) for s in input_shape]]

    # Reshape and sum over appropriate dimensions
    X = X.reshape(sum(zip(final_shape, D), ()))
    return agg(X, tuple([i for i in range(X.ndim) if i % 2 == 1]))


def regrid2(X, D, agg=sum):
    """Regrid an N-dimensional matrix using iterators

    Parameters
    ----------
    X : array_like
        N-dimensional data to regrid
    D : int or array-like
        The step size of the regridding.  D must be either an integer
        or an array of integers of length X.ndim
    agg : function (optional)
        aggregate operation to use (default=sum). Must take an interable as an
        argument, and return the element-wise aggregate of numpy arrays returned
        by the iterator.

    Returns
    -------
    X_agg : numpy array
        for X.shape = (N1, N2, N3...) and D = (D1, D2, D3...) then the shape
        of X_agg is (N1 // D1, N2 // D2, N3 // D3...)

    Note
    ----
    If Ni is not a multiple of Di, then it will be truncated from the right.
    """
    X = np.asarray(X)
    D = np.zeros(X.ndim, dtype=int) + np.asarray(D).astype(int)
    N = np.asarray(X.shape)
    final_shape = N // D
    input_shape = D * (N // D)

    indices = itertools.product(*(xrange(Di) for Di in D))
    slices = [[slice(i, i + shape_i, Di)
               for i, shape_i, Di in itertools.izip(ind, input_shape, D)]
              for ind in indices]

    return agg(X[slc] for slc in slices)


if __name__ == '__main__':
    X = np.random.random((211, 331, 161))
    V1 = regrid(X, (4, 6, 5))
    V2 = regrid2(X, (4, 6, 5))

    print np.all(V1 == V2)
