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
