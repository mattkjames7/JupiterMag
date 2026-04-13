import JupiterMag as jm
import numpy as np
import pytest


def _assert_cfg_matches(actual, expected):
    assert actual.keys() == expected.keys()

    for key, expected_value in expected.items():
        actual_value = actual[key]
        if isinstance(expected_value, (float, np.floating)):
            assert actual_value == pytest.approx(expected_value, rel=0.0, abs=1e-12)
        elif isinstance(expected_value, (bool, np.bool_)):
            assert bool(actual_value) is bool(expected_value)
        else:
            assert actual_value == expected_value


def test_con2020_config_unknown_keyword_warns_and_is_ignored(capsys):
    default_cfg = dict(jm.Con2020.Config("default"))

    try:
        actual_cfg = jm.Con2020.Config(not_a_real_parameter=123)
        captured = capsys.readouterr()

        assert "Keyword argument not_a_real_parameter unrecognized, ignoring." in captured.out
        _assert_cfg_matches(actual_cfg, default_cfg)
    finally:
        jm.Con2020.Config("default")
