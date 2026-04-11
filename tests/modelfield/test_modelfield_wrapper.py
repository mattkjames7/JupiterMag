import JupiterMag as jm
import numpy as np
import numpy.testing as npt
import pytest


FIELD_POSITIONS = (
    np.array([10.0, 15.0, 20.0]),
    np.array([0.0, -5.0, 3.0]),
    np.array([1.0, 2.0, -4.0]),
)


@pytest.mark.parametrize("internal_model", ["jrm09", "vip4"])
def test_modelfield_internal_only_matches_internal_field(internal_model):
    default_cfg = jm.Internal.Config("default", Model=internal_model)

    try:
        jm.Internal.Config(**default_cfg)
        expected = jm.Internal.Field(*FIELD_POSITIONS)

        actual = jm.ModelField(*FIELD_POSITIONS, IntModel=internal_model, ExtModel=None)

        for i in range(3):
            npt.assert_allclose(actual[i], expected[i], rtol=1e-10, atol=1e-10)
    finally:
        jm.Internal.Config("default", Model=internal_model)


def test_modelfield_external_only_matches_con2020_field():
    default_cfg = jm.Con2020.Config("default")

    try:
        jm.Con2020.Config(**default_cfg)
        expected = jm.Con2020.Field(*FIELD_POSITIONS)

        actual = jm.ModelField(*FIELD_POSITIONS, IntModel=None, ExtModel="con2020")

        for i in range(3):
            npt.assert_allclose(actual[i], expected[i], rtol=1e-10, atol=1e-10)
    finally:
        jm.Con2020.Config("default")


def test_modelfield_scalar_inputs_return_length_one_arrays():
    default_cfg = jm.Internal.Config("default", Model="jrm09")

    try:
        jm.Internal.Config(**default_cfg)
        result = jm.ModelField(10.0, -3.0, 5.0, IntModel="jrm09", ExtModel=None)

        assert len(result) == 3
        for component in result:
            assert isinstance(component, np.ndarray)
            assert component.shape == (1,)
    finally:
        jm.Internal.Config("default", Model="jrm09")
