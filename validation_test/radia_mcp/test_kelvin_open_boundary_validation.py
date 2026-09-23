from __future__ import annotations

import pytest

from radia_mcp.radia_ngsolve.fem_bem_coupling import kelvin_twosphere_shell_dipole
from radia_mcp.radia_ngsolve.kelvin_open_boundary_validation import (
    run_kelvin_open_boundary_validation,
)


@pytest.fixture(scope="module")
def verified_result():
    return run_kelvin_open_boundary_validation({})


def test_genuine_two_sphere_periodic_kelvin_solve_is_digest_bound(verified_result):
    assert verified_result["status"] == "verified"
    assert verified_result["pass"] is True
    assert max(row["rel_err"] for row in verified_result["mode_rows_2d"]) < 0.005
    three = verified_result["three_dimensional"]
    assert three["rel_err"] < 0.025
    assert three["free_residual_inf"] < 1.0e-10
    assert len(three["mesh_sha256"]) == len(three["operator_sha256"]) == 64
    assert verified_result["validated_capabilities"] == [
        "kelvin_open_boundary",
        "open_boundary",
    ]


def test_periodic_hash_direction_regression_and_curve_orders():
    result = kelvin_twosphere_shell_dipole(maxh=0.45, order=1, curve_order=2)
    assert result["mesh"]["periodic_map"] == "kelvin_ext_to_kelvin_int_inverse_translation"
    assert result["mesh"]["curve_order"] == 2
    assert result["free_residual_inf"] < 1.0e-10


@pytest.mark.parametrize(
    "case, message",
    [
        ({"method": "pml"}, "kelvin_transform"),
        ({"pml": True}, "PML"),
        ({"physics_regime": "time_harmonic_wave"}, "static-only"),
        ({"wave_boundary_inference": "reuse_for_waves"}, "forbidden"),
        ({"inner_radius": 1.0, "outer_radius": 1.0}, "smaller"),
        ({"offset": 1.5}, "disjoint"),
    ],
)
def test_invalid_static_policy_or_geometry_fails_closed(case, message):
    with pytest.raises(ValueError, match=message):
        run_kelvin_open_boundary_validation(case)


def test_stale_expected_result_identity_is_rejected(verified_result):
    stale = run_kelvin_open_boundary_validation({"expected_result_sha256": "0" * 64})
    assert stale["pass"] is False
    assert stale["issues"] == ["expected_result_identity_matches"]
    assert len(stale["result_sha256"]) == 64
