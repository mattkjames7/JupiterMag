import JupiterMag as jm
import numpy as np
import numpy.testing as npt
import pytest


@pytest.mark.parametrize(
    "x, y, z, xt, xp",
    [
        (10.0, -3.0, 5.0, 9.3, 155.8),
        (
            np.array([1.0, 5.0, -2.0, 7.5]),
            np.array([0.0, -4.0, 3.0, 1.5]),
            np.array([2.0, 6.0, -8.0, 0.5]),
            9.3,
            155.8,
        ),
        (
            np.array([3.0, -1.0, 2.5]),
            np.array([4.0, 0.5, -6.0]),
            np.array([1.0, 7.0, -2.0]),
            25.0,
            45.0,
        ),
    ],
)
def test_coordconv_roundtrip_siii_mag_siii(x, y, z, xt, xp):
    mx, my, mz = jm.CoordConv.SIIItoMag(x, y, z, xt, xp)
    rx, ry, rz = jm.CoordConv.MagtoSIII(mx, my, mz, xt, xp)

    npt.assert_allclose(np.asarray(rx), np.asarray(x), rtol=1e-12, atol=1e-12)
    npt.assert_allclose(np.asarray(ry), np.asarray(y), rtol=1e-12, atol=1e-12)
    npt.assert_allclose(np.asarray(rz), np.asarray(z), rtol=1e-12, atol=1e-12)


@pytest.mark.parametrize(
    "coords",
    [
        (10.0, -3.0, 5.0),
        (
            np.array([1.0, -2.0, 3.0]),
            np.array([4.0, 5.0, -6.0]),
            np.array([7.0, -8.0, 9.0]),
        ),
    ],
)
def test_coordconv_zero_tilt_is_identity(coords):
    x, y, z = coords

    mx, my, mz = jm.CoordConv.SIIItoMag(x, y, z, 0.0, 0.0)
    sx, sy, sz = jm.CoordConv.MagtoSIII(x, y, z, 0.0, 0.0)

    npt.assert_allclose(np.asarray(mx), np.asarray(x), rtol=0.0, atol=1e-12)
    npt.assert_allclose(np.asarray(my), np.asarray(y), rtol=0.0, atol=1e-12)
    npt.assert_allclose(np.asarray(mz), np.asarray(z), rtol=0.0, atol=1e-12)
    npt.assert_allclose(np.asarray(sx), np.asarray(x), rtol=0.0, atol=1e-12)
    npt.assert_allclose(np.asarray(sy), np.asarray(y), rtol=0.0, atol=1e-12)
    npt.assert_allclose(np.asarray(sz), np.asarray(z), rtol=0.0, atol=1e-12)
