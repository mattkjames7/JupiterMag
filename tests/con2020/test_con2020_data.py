import json
import sys
from pathlib import Path

import numpy as np
import numpy.testing as npt
import pytest

CURRENT_DIR = Path(__file__).resolve().parent
sys.path.append(str(CURRENT_DIR.parent))
from common import get_con2020_field  # noqa: E402

DATA_FILE = Path(__file__).resolve().parents[1] / "data" / "con2020_field_data.json"


def _load_test_data():
    with DATA_FILE.open("r", encoding="utf-8") as f:
        return json.load(f)


TEST_CASES = _load_test_data()


def _case_id(case):
    cfg = case["input"]["args"][3]
    return (
        f"Eq={cfg['equation_type']}-"
        f"Cin={cfg['CartesianIn']}-"
        f"Cout={cfg['CartesianOut']}-"
        f"Edw={cfg['Edwards']}-"
        f"Az={cfg['azfunc']}-"
        f"Sm={cfg['Smooth']}"
    )


@pytest.mark.parametrize("case", TEST_CASES, ids=_case_id)
def test_con2020_saved_data(case):
    args = case["input"]["args"]
    expected = [np.asarray(component) for component in case["output"]["result"]]

    result = get_con2020_field(*args)
    actual = [np.asarray(component) for component in result]

    assert len(actual) == 3
    for i in range(3):
        npt.assert_allclose(actual[i], expected[i], rtol=1e-10, atol=1e-10)
