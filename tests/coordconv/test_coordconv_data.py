import json
from pathlib import Path

import numpy as np
import numpy.testing as npt
import pytest

import JupiterMag as jm

DATA_FILE = Path(__file__).resolve().parents[1] / "data" / "coordconv_data.json"


def _load_test_data():
    with DATA_FILE.open("r", encoding="utf-8") as f:
        return json.load(f)


TEST_CASES = _load_test_data()

FUNCTIONS = {
    "MagtoSIII": jm.CoordConv.MagtoSIII,
    "SIIItoMag": jm.CoordConv.SIIItoMag,
}


def _case_id(case):
    fn = case["function"]
    x, y, z, xt, xp = case["input"]["args"]
    return f"{fn}-xt={xt}-xp={xp}-x={x}-y={y}-z={z}"


@pytest.mark.parametrize("case", TEST_CASES, ids=_case_id)
def test_coordconv_saved_data(case):
    function_name = case["function"]
    args = case["input"]["args"]
    kwargs = case["input"]["kwargs"]
    expected = np.asarray(case["output"]["result"])

    result = FUNCTIONS[function_name](*args, **kwargs)
    actual = np.asarray(result)

    npt.assert_allclose(actual, expected, rtol=1e-12, atol=1e-12)
