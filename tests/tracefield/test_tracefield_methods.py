import JupiterMag as jm
import numpy as np


def test_tracefield_wrapextfuncs_accepts_string():
    trace = jm.TraceField.__new__(jm.TraceField)

    n_ext, ext_ptr = trace._WrapExtFuncs("Con2020")

    assert int(n_ext) == 1
    assert ext_ptr[0] == b"Con2020"


def test_tracefield_wrapextfuncs_accepts_list():
    trace = jm.TraceField.__new__(jm.TraceField)

    n_ext, ext_ptr = trace._WrapExtFuncs(["Con2020", "none"])

    assert int(n_ext) == 2
    assert list(ext_ptr[:2]) == [b"Con2020", b"none"]


def test_tracefield_storetime_with_date_and_ut_sets_time_fields():
    trace = jm.TraceField.__new__(jm.TraceField)
    date = np.array([20200101, 20200101], dtype="int32")
    ut = np.array([0.0, 1.5], dtype="float64")

    trace.StoreTime(Date=date, ut=ut)

    assert trace.Time is True
    np.testing.assert_array_equal(trace.Date, date)
    np.testing.assert_array_equal(trace.ut, ut)
    assert trace.utc.shape == ut.shape


def test_tracefield_storetime_none_leaves_time_disabled():
    trace = jm.TraceField.__new__(jm.TraceField)

    trace.StoreTime(Time=None)

    assert trace.Time is False
