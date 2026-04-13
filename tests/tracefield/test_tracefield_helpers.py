import JupiterMag as jm
import numpy as np
import numpy.testing as npt


def _make_trace_dict():
    n = 2
    max_len = 4
    nalpha = 2

    x = np.array([[1.0, 2.0, 3.0, np.nan], [4.0, 5.0, np.nan, np.nan]])
    y = x + 10.0
    z = x + 20.0
    bx = x + 30.0
    by = x + 40.0
    bz = x + 50.0
    s = x + 60.0
    r = x + 70.0
    rnorm = x + 80.0
    halpha = np.arange(n * nalpha * max_len, dtype="float64").reshape(n, nalpha, max_len)

    return {
        "n": n,
        "MaxLen": max_len,
        "nalpha": np.int32(nalpha),
        "alpha": np.array([0.0, 30.0]),
        "nstep": np.array([3, 2], dtype="int32"),
        "x": x,
        "y": y,
        "z": z,
        "Bx": bx,
        "By": by,
        "Bz": bz,
        "s": s,
        "R": r,
        "Rnorm": rnorm,
        "halpha": halpha,
        "IntModelCode": b"jrm09",
        "ExtModelCode": [b"con2020"],
        "traceRegion": np.zeros((n, max_len), dtype="int32"),
        "FP": np.zeros((n, 49), dtype="float64"),
    }


def _make_trace_with_footprints():
    data = _make_trace_dict()
    data["FP"] = np.arange(data["n"] * 49, dtype="float64").reshape(data["n"], 49)
    return jm.TraceField(data)


def test_tracefield_init_from_dict_preserves_data():
    data = _make_trace_dict()
    trace = jm.TraceField(data)

    assert trace.n == data["n"]
    assert trace.MaxLen == data["MaxLen"]
    npt.assert_array_equal(trace.nstep, data["nstep"])
    npt.assert_array_equal(trace.x, data["x"])
    npt.assert_array_equal(trace.halpha, data["halpha"])


def test_tracefield_tracedict_removes_nan_padding_and_ctypes_fields():
    trace = jm.TraceField(_make_trace_dict())

    out = trace.TraceDict(RemoveNAN=True)

    assert "IntModelCode" not in out
    assert "ExtModelCode" not in out

    assert out["x"].dtype == object
    assert out["halpha"].dtype == object

    npt.assert_array_equal(out["x"][0], np.array([1.0, 2.0, 3.0]))
    npt.assert_array_equal(out["x"][1], np.array([4.0, 5.0]))
    npt.assert_array_equal(out["Rnorm"][0], np.array([81.0, 82.0, 83.0]))
    npt.assert_array_equal(out["halpha"][0, 0], np.array([0.0, 1.0, 2.0]))
    npt.assert_array_equal(out["halpha"][0, 1], np.array([4.0, 5.0, 6.0]))
    npt.assert_array_equal(out["halpha"][1, 0], np.array([8.0, 9.0]))


def test_tracefield_unpackfootprints_maps_fp_columns_to_named_fields():
    trace = _make_trace_with_footprints()

    trace._UnpackFootprints()

    assert trace.ionosphere.xn3.shape == (trace.n,)
    assert trace.surface.lonn.shape == (trace.n,)
    assert trace.equator.fllen.shape == (trace.n,)

    npt.assert_array_equal(trace.ionosphere.xn3, trace.FP[:, 0])
    npt.assert_array_equal(trace.ionosphere.latn, trace.FP[:, 31])
    npt.assert_array_equal(trace.surface.xn3, trace.FP[:, 12])
    npt.assert_array_equal(trace.surface.mlats, trace.FP[:, 45])
    npt.assert_array_equal(trace.equator.x3, trace.FP[:, 24])
    npt.assert_array_equal(trace.equator.fllen, trace.FP[:, 48])


def test_tracefield_gettrace_returns_trimmed_trace_arrays():
    trace = jm.TraceField(_make_trace_dict())

    x, y, z, bx, by, bz, r, rnorm, s, h = trace.GetTrace(0)

    npt.assert_array_equal(x, np.array([1.0, 2.0, 3.0]))
    npt.assert_array_equal(y, np.array([11.0, 12.0, 13.0]))
    npt.assert_array_equal(z, np.array([21.0, 22.0, 23.0]))
    npt.assert_array_equal(bx, np.array([31.0, 32.0, 33.0]))
    npt.assert_array_equal(by, np.array([41.0, 42.0, 43.0]))
    npt.assert_array_equal(bz, np.array([51.0, 52.0, 53.0]))
    npt.assert_array_equal(r, np.array([71.0, 72.0, 73.0]))
    npt.assert_array_equal(rnorm, np.array([81.0, 82.0, 83.0]))
    npt.assert_array_equal(s, np.array([61.0, 62.0, 63.0]))
    assert h.shape == (trace.nalpha, trace.nstep[0])
    npt.assert_array_equal(h, np.array([[0.0, 1.0, 2.0], [4.0, 5.0, 6.0]]))
