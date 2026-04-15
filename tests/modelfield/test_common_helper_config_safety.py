import importlib
from pathlib import Path
import sys

import JupiterMag as jm
import numpy as np

CURRENT_DIR = Path(__file__).resolve().parent
COMMON_DIR = CURRENT_DIR.parent
if str(COMMON_DIR) not in sys.path:
    sys.path.append(str(COMMON_DIR))

import common  # noqa: E402


def test_common_import_does_not_mutate_global_configs():
    internal_before = dict(jm.Internal.Config())
    con2020_before = dict(jm.Con2020.Config())

    importlib.reload(common)

    internal_after = dict(jm.Internal.Config())
    con2020_after = dict(jm.Con2020.Config())

    assert internal_after == internal_before
    assert con2020_after == con2020_before


def test_get_model_field_restores_global_configs():
    internal_before = dict(jm.Internal.Config())
    con2020_before = dict(jm.Con2020.Config())

    x = np.array([10.0, 12.0], dtype="float64")
    y = np.array([1.0, -2.0], dtype="float64")
    z = np.array([3.0, 4.0], dtype="float64")

    common.get_model_field(x, y, z, IntModel="jrm09", ExtModel="con2020")

    internal_after = dict(jm.Internal.Config())
    con2020_after = dict(jm.Con2020.Config())

    assert internal_after == internal_before
    assert con2020_after == con2020_before


def test_get_internal_field_restores_global_config():
    internal_before = dict(jm.Internal.Config())

    common.get_internal_field(
        np.array([10.0, 12.0], dtype="float64"),
        np.array([1.0, -2.0], dtype="float64"),
        np.array([3.0, 4.0], dtype="float64"),
        model="vip4",
    )

    internal_after = dict(jm.Internal.Config())
    assert internal_after == internal_before


def test_get_con2020_field_restores_global_config():
    con2020_before = dict(jm.Con2020.Config())

    cfg = dict(jm.Con2020.Config("default"))
    common.get_con2020_field(
        np.array([10.0, 12.0], dtype="float64"),
        np.array([1.0, -2.0], dtype="float64"),
        np.array([3.0, 4.0], dtype="float64"),
        cfg,
    )

    con2020_after = dict(jm.Con2020.Config())
    assert con2020_after == con2020_before
