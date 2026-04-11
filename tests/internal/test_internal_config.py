import JupiterMag as jm
import numpy as np
import numpy.testing as npt
import pytest


FIELD_POSITIONS = (
    np.array([0.0, 5.0, 10.0, 20.0, 10.0, 10.0, 10.0, 10.0, 20.0, 20.0]),
    np.array([5.0, 0.0, 0.0, 0.0, 5.0, 10.0, 20.0, 10.0, 5.0, 10.0]),
    np.array([20.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 5.0, 10.0, 20.0]),
)

MODEL_DEGREES = {
    "jrm09": {"full": 10, "reduced": 9},
    "vip4": {"full": 4, "reduced": 3},
}

DEFAULT_INTERNAL_CFG = {
    "Model": "jrm09",
    "CartesianIn": True,
    "CartesianOut": True,
    "Degree": 0,
}


def _default_config(model="jrm09"):
    return jm.Internal.Config("default", Model=model)


def _get_internal_field(cfg):
    jm.Internal.Config(**cfg)
    return jm.Internal.Field(*FIELD_POSITIONS)


def _get_modified_value(parameter, value):
    if isinstance(value, (bool, np.bool_)):
        return not value

    if parameter == "Degree":
        return max(1, int(value) - 1) if int(value) > 1 else 1

    if parameter == "Model":
        return "vip4" if str(value).lower() != "vip4" else "jrm09"

    raise ValueError(f"No alternate value defined for {parameter}")


def _assert_cfg_matches(actual, expected):
    assert actual.keys() == expected.keys()

    for key, expected_value in expected.items():
        actual_value = actual[key]
        if isinstance(expected_value, (bool, np.bool_)):
            assert bool(actual_value) is bool(expected_value)
        elif isinstance(expected_value, (int, np.integer)):
            assert int(actual_value) == int(expected_value)
        else:
            assert actual_value == expected_value


@pytest.mark.parametrize("model", sorted(MODEL_DEGREES))
def test_internal_config_degree_change_modifies_field_and_restore_default(model):
    model_degrees = MODEL_DEGREES[model]
    baseline_cfg = {
        "Model": model,
        "CartesianIn": True,
        "CartesianOut": True,
        "Degree": model_degrees["full"],
    }

    try:
        original = _get_internal_field(baseline_cfg)

        modified_cfg = dict(baseline_cfg)
        modified_cfg["Degree"] = model_degrees["reduced"]
        updated = _get_internal_field(modified_cfg)

        original_vectors = np.column_stack(original)
        updated_vectors = np.column_stack(updated)
        vector_delta = np.linalg.norm(updated_vectors - original_vectors, axis=1)

        assert np.any(vector_delta > 0.0)
        assert not np.allclose(updated_vectors, original_vectors, rtol=1e-12, atol=1e-12)

        restored = _get_internal_field(baseline_cfg)
        for i in range(3):
            npt.assert_allclose(restored[i], original[i], rtol=1e-10, atol=1e-10)
    finally:
        jm.Internal.Config("default", Model=model)


@pytest.mark.parametrize("model", sorted(MODEL_DEGREES))
@pytest.mark.parametrize("parameter", sorted(DEFAULT_INTERNAL_CFG))
def test_internal_config_changes_persist(model, parameter):
    default_cfg = dict(_default_config(model))

    try:
        current_cfg = jm.Internal.Config("default", Model=model)
        _assert_cfg_matches(current_cfg, default_cfg)

        updated_value = _get_modified_value(parameter, default_cfg[parameter])

        if parameter == "Model":
            expected_cfg = dict(_default_config(updated_value))
        else:
            expected_cfg = dict(default_cfg)
            expected_cfg[parameter] = updated_value

        updated_cfg = jm.Internal.Config(**{parameter: updated_value})
        _assert_cfg_matches(updated_cfg, expected_cfg)
    finally:
        jm.Internal.Config("default", Model=model)