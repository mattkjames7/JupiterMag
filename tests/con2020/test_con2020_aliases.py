import JupiterMag as jm
import numpy as np
import pytest

LONGNAME_ALIASES = {
    "mu_i_div2__current_parameter_nT": "mu_i",
    "i_rho__radial_current_MA": "i_rho",
    "r0__inner_rj": "r0",
    "r1__outer_rj": "r1",
    "d__cs_half_thickness_rj": "d",
    "xt__cs_tilt_degs": "xt",
    "xp__cs_rhs_azimuthal_angle_of_tilt_degs": "xp",
}


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


@pytest.mark.parametrize("long_name, short_name", sorted(LONGNAME_ALIASES.items()))
def test_con2020_config_long_name_alias_persists(long_name, short_name):
    default_cfg = dict(jm.Con2020.Config("default"))

    try:
        updated_value = default_cfg[short_name] + 1.25
        expected_cfg = dict(default_cfg)
        expected_cfg[short_name] = updated_value

        updated_cfg = jm.Con2020.Config(**{long_name: updated_value})
        _assert_cfg_matches(updated_cfg, expected_cfg)
    finally:
        jm.Con2020.Config("default")


@pytest.mark.parametrize("long_name, short_name", sorted(LONGNAME_ALIASES.items()))
def test_con2020_config_long_name_alias_matches_short_name(long_name, short_name):
    default_cfg = dict(jm.Con2020.Config("default"))
    updated_value = default_cfg[short_name] + 2.5

    try:
        cfg_from_long_name = jm.Con2020.Config(**{long_name: updated_value})
        jm.Con2020.Config("default")
        cfg_from_short_name = jm.Con2020.Config(**{short_name: updated_value})

        _assert_cfg_matches(cfg_from_long_name, cfg_from_short_name)
    finally:
        jm.Con2020.Config("default")
