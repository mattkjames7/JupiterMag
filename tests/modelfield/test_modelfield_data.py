import json
from pathlib import Path

import numpy as np
import numpy.testing as npt
import pytest

import sys
CURRENT_DIR = Path(__file__).resolve().parent
sys.path.append(str(CURRENT_DIR.parent))
from common import get_model_field


DATA_FILE = Path(__file__).resolve().parents[1] / "data" / "modelfield_data.json"


def _load_test_data():
    with DATA_FILE.open("r", encoding="utf-8") as f:
        return json.load(f)


TEST_CASES = _load_test_data()


def _case_id(case):
    kwargs = case["input"]["kwargs"]
    int_model = kwargs["IntModel"]
    ext_model = kwargs["ExtModel"]
    cart_in = kwargs["CartIn"]
    cart_out = kwargs["CartOut"]
    return f"Int={int_model}-Ext={ext_model}-CartIn={cart_in}-CartOut={cart_out}"


@pytest.mark.parametrize("case", TEST_CASES, ids=_case_id)
def test_modelfield_saved_data(case):
    args = case["input"]["args"]
    kwargs = case["input"]["kwargs"]
    expected = [np.asarray(component) for component in case["output"]["result"]]

    result = get_model_field(*args, **kwargs)
    actual = [np.asarray(component) for component in result]

    assert len(actual) == 3
    for i in range(3):
        npt.assert_allclose(actual[i], expected[i], rtol=1e-10, atol=1e-10)
