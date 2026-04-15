import json
import sys
from pathlib import Path

import numpy as np
import numpy.testing as npt
import pytest

CURRENT_DIR = Path(__file__).resolve().parent
sys.path.append(str(CURRENT_DIR.parent))
from common import get_trace_footprints  # noqa: E402

DATA_FILE = Path(__file__).resolve().parents[1] / "data" / "trace_data.json"


def _load_test_data():
    with DATA_FILE.open("r", encoding="utf-8") as f:
        return json.load(f)


TEST_CASES = _load_test_data()


def _case_id(case):
    x, y, z, int_model, ext_model, con2020_cfg = case["input"]["args"]
    eq_type = con2020_cfg["equation_type"]
    r = float(np.sqrt(x * x + y * y + z * z))
    return f"Int={int_model}-Ext={ext_model}-Eq={eq_type}-r={r:.1f}"


def _assert_nested_allclose(actual, expected):
    if isinstance(expected, dict):
        assert isinstance(actual, dict)
        assert set(actual.keys()) == set(expected.keys())
        for key in expected:
            _assert_nested_allclose(actual[key], expected[key])
        return

    actual_array = np.asarray(actual)
    expected_array = np.asarray(expected)
    npt.assert_allclose(actual_array, expected_array, rtol=1e-10, atol=1e-10)


@pytest.mark.parametrize("case", TEST_CASES, ids=_case_id)
def test_tracefield_saved_data(case):
    args = case["input"]["args"]
    kwargs = case["input"]["kwargs"]
    expected = case["output"]["result"]

    actual = get_trace_footprints(*args, **kwargs)

    _assert_nested_allclose(actual, expected)
