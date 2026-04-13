import JupiterMag as jm
import numpy as np
import numpy.testing as npt


FIELD_POSITIONS = (
    np.array([0.0, 5.0, 10.0, 20.0]),
    np.array([5.0, 0.0, 0.0, 0.0]),
    np.array([20.0, 0.0, 0.0, 0.0]),
)


def test_internal_field_maxdeg_override_does_not_modify_stored_config():
    try:
        stored_cfg = jm.Internal.Config("default", Model="jrm09", Degree=10)
        expected_default = jm.Internal.Field(*FIELD_POSITIONS)
        expected_explicit = jm.Internal.Field(*FIELD_POSITIONS, MaxDeg=10)
        reduced = jm.Internal.Field(*FIELD_POSITIONS, MaxDeg=5)

        for i in range(3):
            npt.assert_allclose(expected_default[i], expected_explicit[i], rtol=1e-10, atol=1e-10)

        assert not np.allclose(np.column_stack(expected_default), np.column_stack(reduced), rtol=1e-12, atol=1e-12)

        cfg_after = jm.Internal.Config()
        assert int(cfg_after["Degree"]) == int(stored_cfg["Degree"])
        assert cfg_after["Model"] == stored_cfg["Model"]
    finally:
        jm.Internal.Config("default", Model="jrm09")


def test_internal_config_unknown_keyword_warns_and_is_ignored(capsys):
    default_cfg = dict(jm.Internal.Config("default", Model="jrm09"))

    try:
        actual_cfg = jm.Internal.Config(not_a_real_parameter=123)
        captured = capsys.readouterr()

        assert "Keyword argument not_a_real_parameter unrecognized, ignoring." in captured.out
        assert actual_cfg["Model"] == default_cfg["Model"]
        assert bool(actual_cfg["CartesianIn"]) is bool(default_cfg["CartesianIn"])
        assert bool(actual_cfg["CartesianOut"]) is bool(default_cfg["CartesianOut"])
        assert int(actual_cfg["Degree"]) == int(default_cfg["Degree"])
    finally:
        jm.Internal.Config("default", Model="jrm09")
