import JupiterMag as jm
import numpy as np
import numpy.testing as npt
import pytest

FIELD_POSITIONS = (
    np.array([10.0, 15.0, -20.0, 30.0]),
    np.array([1.0, -3.0, 5.0, 7.0]),
    np.array([2.0, 4.0, -6.0, 8.0]),
)


def _default_config():
    return jm.Con2020.Config("default")


def _get_con2020_field(cfg):
    jm.Con2020.Config(**cfg)
    return jm.Con2020.Field(*FIELD_POSITIONS)


def _get_modified_value(key, value):
    if isinstance(value, (bool, np.bool_)):
        return not value

    if isinstance(value, str):
        if key == "equation_type":
            return "analytic" if value != "analytic" else "integral"
        if key == "azfunc":
            return "lmic" if value != "lmic" else "connerney"
        raise ValueError(f"No alternate value defined for {key}")

    return value + 1.25


def _assert_cfg_matches(actual, expected):
    assert actual.keys() == expected.keys()

    for key, expected_value in expected.items():
        actual_value = actual[key]
        if isinstance(expected_value, (float, np.floating)):
            assert actual_value == pytest.approx(expected_value, rel=0.0, abs=1e-12)
        else:
            assert actual_value == expected_value


def test_con2020_config_changes_modify_field_and_restore_default():
    default_cfg = dict(_default_config())

    try:
        original = _get_con2020_field(default_cfg)

        modified_cfg = dict(default_cfg)
        modified_cfg["i_rho"] = default_cfg["i_rho"] + 5.0
        modified_cfg["mu_i"] = default_cfg["mu_i"] + 25.0
        updated = _get_con2020_field(modified_cfg)

        original_vectors = np.column_stack(original)
        updated_vectors = np.column_stack(updated)
        vector_delta = np.linalg.norm(updated_vectors - original_vectors, axis=1)

        assert np.all(vector_delta > 0.0)

        restored = _get_con2020_field(default_cfg)
        for i in range(3):
            npt.assert_allclose(restored[i], original[i], rtol=1e-10, atol=1e-10)
    finally:
        jm.Con2020.Config("default")


@pytest.mark.parametrize("parameter", sorted(_default_config()))
def test_con2020_config_changes_persist(parameter):
    default_cfg = dict(_default_config())

    try:
        current_cfg = jm.Con2020.Config("default")
        _assert_cfg_matches(current_cfg, default_cfg)

        updated_value = _get_modified_value(parameter, default_cfg[parameter])
        expected_cfg = dict(default_cfg)
        expected_cfg[parameter] = updated_value

        updated_cfg = jm.Con2020.Config(**{parameter: updated_value})
        _assert_cfg_matches(updated_cfg, expected_cfg)
    finally:
        jm.Con2020.Config("default")
