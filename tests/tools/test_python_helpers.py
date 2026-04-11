import numpy as np

from JupiterMag.Tools.GetLegendHandLab import GetLegendHandLab
from JupiterMag.Tools.JupiterOval import JupiterOvalNorth, JupiterOvalSouth
from JupiterMag._ptr2D import _ptr2D


class _FakeText:
    def __init__(self, text):
        self._text = text

    def get_text(self):
        return self._text


class _FakeLegend:
    def __init__(self, handles, labels):
        self.legendHandles = handles
        self.texts = [_FakeText(label) for label in labels]


class _FakeAxes:
    def __init__(self, legend=None):
        self.legend_ = legend


def test_getlegendhandlab_returns_empty_lists_without_legend():
    handles, labels = GetLegendHandLab(_FakeAxes())

    assert handles == []
    assert labels == []


def test_getlegendhandlab_returns_handles_and_labels_from_legend():
    handles_in = [object(), object()]
    labels_in = ["north", "south"]

    handles, labels = GetLegendHandLab(_FakeAxes(_FakeLegend(handles_in, labels_in)))

    assert handles == handles_in
    assert labels == labels_in


def test_ptr2d_returns_row_start_addresses():
    arr = np.arange(12, dtype="float64").reshape(3, 4)

    ptrs = _ptr2D(arr)

    assert ptrs.shape == (3,)
    assert ptrs.dtype == np.uintp
    np.testing.assert_array_equal(np.diff(ptrs), np.full(2, arr.strides[0], dtype=np.uintp))


def test_jupiter_oval_outputs_are_finite_and_have_expected_signs():
    north_lon, north_lat = JupiterOvalNorth()
    south_lon, south_lat = JupiterOvalSouth()

    assert north_lon.shape == north_lat.shape
    assert south_lon.shape == south_lat.shape
    assert np.all(np.isfinite(north_lon))
    assert np.all(np.isfinite(north_lat))
    assert np.all(np.isfinite(south_lon))
    assert np.all(np.isfinite(south_lat))
    assert np.all(north_lat > 0.0)
    assert np.all(south_lat < 0.0)
