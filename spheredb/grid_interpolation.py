import numpy as np
import matplotlib.pyplot as plt

__all__ = ['GridInterpolation']


class GridInterpolation(object):
    """Interpolation on a 2D regular grid

    Parameters
    ----------
    z : array-like
        2-dimensional array
    origin, step : array-like, shape = z.ndim
        The grid specification.
        x = origin[0] + step[0] * i
        y = origin[1] + step[1] * j
        f(x, y) = z

    Calling
    -------
    Call with an input array X, of shape (n_samples, 2)

    Examples
    --------
    [TODO]
    """
    def __init__(self, z, origin, step):
        z, origin, step = map(np.asarray, (z, origin, step))
        if origin.shape != (z.ndim,) or step.shape != (z.ndim,):
            raise ValueError("Number of dimensions in z should match "
                             "number of items in origin and step")

        if z.ndim != 2:
            raise NotImplementedError("only implemented for 2D case")

        self.z = z
        self.origin = origin.astype(float)
        self.step = step.astype(float)
        self.volume = np.prod(step)
        self.Nsteps = z.shape

    def __call__(self, X):
        X = np.asarray(X)
        output_shape = X.shape[:-1]
        X = np.atleast_2d(X)

        if X.shape[-1:] != self.origin.shape:
            raise ValueError("dimension of x must match dimension of input")

        ind = (X - self.origin) / self.step
        ind_floor = np.floor(ind).astype(int)
        ind_ceil = np.ceil(ind).astype(int)

        results = np.empty(X.shape[:-1])
        results.fill(np.nan)

        in_bounds = np.logical_and.reduce((ind_floor >= 0)
                                          & (ind_ceil < self.z.shape), -1)

        X = X[in_bounds]
        ind_floor = ind_floor[in_bounds]
        ind_ceil = ind_ceil[in_bounds]

        f11 = self.z[ind_floor[:, 0], ind_floor[:, 1]]
        f12 = self.z[ind_floor[:, 0], ind_ceil[:, 1]]
        f21 = self.z[ind_ceil[:, 0], ind_floor[:, 1]]
        f22 = self.z[ind_ceil[:, 0], ind_ceil[:, 1]]

        X1 = (X - self.origin) % self.step
        X0 = self.step - X1

        X0 = X0.T
        X1 = X1.T

        results[in_bounds] = (f11 * X0[0] * X0[1] +
                              f12 * X0[0] * X1[1] +
                              f21 * X1[0] * X0[1] +
                              f22 * X1[0] * X1[1]) / self.volume

        return results.reshape(output_shape)


def plot_test():
    import matplotlib.pyplot as plt

    x = np.linspace(0, 10, 20)
    y = np.linspace(0, 10, 19)
    dx = x[1] - x[0]
    dy = y[1] - y[0]
    z = ((np.sin(x[:, None]) + np.cos(y - 5)) /
         (1 + np.sqrt((x[:, None] - 2) ** 2 + (y - 3) ** 2)))
    
    I = GridInterpolation(z, [x[0], y[0]], [dx, dy])

    x2 = np.linspace(-1, 11, 40)
    y2 = np.linspace(-1, 11, 35)
    X = np.asarray(np.meshgrid(x2, y2)).transpose((2, 1, 0))
    z2 = I(X.reshape(-1, 2)).reshape(X.shape[:-1])

    plt.subplot(211)
    plt.contourf(x, y, z.T, 50)
    plt.colorbar()
    plt.xlim(0, 10)
    plt.ylim(0, 10)
    
    plt.subplot(212)
    plt.contourf(x2, y2, z2.T, 50)
    plt.colorbar()
    plt.xlim(0, 10)
    plt.ylim(0, 10)
    
    plt.show()


if __name__ == '__main__':
    plot_test()
