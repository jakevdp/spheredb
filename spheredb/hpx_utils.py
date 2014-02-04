"""HEALPix Utilities"""
import numpy as np


def RAdec_to_HPX(RA, dec):
    """Convert RA/dec to healpix

    Parameters
    ----------
    RA, dec : degrees

    Returns
    -------
    x, y : degrees

    See Section 6 of Calabretta & Roukema, Mapping on the HEALPix grid
    """
    H, K = 4.0, 3.0

    RA = np.asarray(RA, dtype=float)
    dec = np.asarray(dec, dtype=float)

    # shift the RA to the range [-180, 180)
    RA = -180 + (RA + 180) % 360

    RA = RA + np.zeros_like(dec)
    dec = dec + np.zeros_like(RA)

    if np.any(RA < -180.) or np.any(RA >= 180.):
        raise ValueError("RA must be in range [-180, 180)")

    if np.any(dec < -90) or np.any(dec > 90):
        raise ValueError("DEC must be in range [-90, 90]")

    x = np.zeros(RA.shape, dtype=float)
    y = np.zeros(dec.shape, dtype=float)

    sindec = np.sin(np.radians(dec))

    dec_cutoff = np.degrees(np.arcsin((K - 1.) / K))
    sigma = np.sqrt(K * (1 - abs(sindec)))
    omega = ((K % 2 > 0) | (dec > 0)).astype(float)
    phi_c = -180. + (180. / H) * (omega + 
                                  2 * np.floor((RA + 180.) * H / 360.
                                               + 0.5 * (1. - omega)))

    upper = (dec > dec_cutoff)
    lower = (dec < -dec_cutoff)
    inner = ~(upper | lower)

    x[upper] = phi_c[upper] + (RA[upper] - phi_c[upper]) * sigma[upper]
    y[upper] = (180. / H) * (0.5 * (K + 1) - sigma[upper])

    x[inner] = RA[inner]
    y[inner] = (K * 90. / H) * sindec[inner]

    x[lower] = phi_c[lower] + (RA[lower] - phi_c[lower]) * sigma[lower]
    y[lower] = -(180. / H) * (0.5 * (K + 1) - sigma[lower])

    return x, y


def HPX_to_RAdec(x, y):
    """Convert RA/dec to healpix

    Parameters
    ----------
    RA, dec : degrees

    Returns
    -------
    x, y : degrees

    See Section 6 of Calabretta & Roukema, Mapping on the HEALPix grid
    """
    H, K = 4.0, 3.0

    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)

    x = x + np.zeros_like(y)
    y = y + np.zeros_like(x)

    RA = np.zeros(x.shape, dtype=float)
    dec = np.zeros(y.shape, dtype=float)

    extreme = (abs(y) >= 90)

    upper = ~extreme & (y > (K - 1.) * 90. / H)
    lower = ~extreme & (y < -(K - 1.) * 90. / H)
    inner = ~(upper | lower | extreme)

    sigma = 0.5 * (K + 1) - abs(y * H) / 180.
    omega = ((K % 2 > 0) | (y > 0)).astype(float)
    x_c = -180. + (2 * np.floor((x + 180.) * H / 360. + 0.5 * (1 - omega))
                   + omega) * 180. / H

    RA[upper] = x_c[upper] + (x[upper] - x_c[upper]) / sigma[upper]
    dec[upper] = np.degrees(np.arcsin(1 - (1. / K) * sigma[upper] ** 2))

    RA[inner] = x[inner]
    dec[inner] = np.degrees(np.arcsin((y[inner] * H) / (90. * K)))

    RA[lower] = x_c[lower] + (x[lower] - x_c[lower]) / sigma[lower]
    dec[lower] = -np.degrees(np.arcsin(1 - (1. / K) * sigma[lower] ** 2))

    RA[extreme] = x[extreme]
    dec[extreme] = y[extreme]

    return RA, dec
