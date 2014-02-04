import numpy as np
from numpy.testing import assert_allclose
from spheredb.hpx_utils import RAdec_to_HPX, HPX_to_RAdec


def test_extreme_RA():
    RA1 = np.linspace(-360, 360, 50)
    RA2 = (RA1 + 180) % 360 - 180
    dec = np.linspace(-90, 90, 25)[1:-1, None]

    x1, y1 = RAdec_to_HPX(RA1, dec)
    x2, y2 = RAdec_to_HPX(RA2, dec)

    assert_allclose(x1, x2)
    assert_allclose(y1, y2)


def test_extreme_dec():
    RA = np.linspace(-180, 180)[:-1]
    dec = np.array([[-90], [90]])

    x, y = RAdec_to_HPX(RA, dec)
    RA2, dec2 = HPX_to_RAdec(x, y)

    # Note that RA is undefined here, and input/output will not match
    assert_allclose(dec + np.zeros_like(RA), dec2)


def test_roundtrip():
    RA = np.linspace(-180., 180., 50)[:-1]
    dec = np.linspace(-90., 90., 25)[1:-1, None]

    x, y = RAdec_to_HPX(RA, dec)
    RA_out, dec_out = HPX_to_RAdec(x, y)

    assert_allclose(RA + np.zeros_like(dec), RA_out)
    assert_allclose(dec + np.zeros_like(RA), dec_out)
