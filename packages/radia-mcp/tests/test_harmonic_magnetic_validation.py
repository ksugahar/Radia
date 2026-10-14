from __future__ import annotations

import pytest

from radia_mcp.radia_ngsolve.harmonic_magnetic_validation import (
    run_harmonic_magnetic_validation,
)


@pytest.fixture(scope="module")
def verified():
    return run_harmonic_magnetic_validation({})


def test_planar_and_axisymmetric_complex_solves_are_verified(verified):
    assert verified["status"] == "verified"
    assert verified["pass"] is True
    assert verified["validated_capabilities"] == [
        "axisymmetric_ac_magnetic",
        "planar_ac_magnetic",
    ]
    assert {row["formulation"] for row in verified["cases"]} == {
        "planar_az",
        "axisymmetric_psi_r_aphi",
    }


def test_errors_residuals_losses_and_identities_are_physical(verified):
    for row in verified["cases"]:
        assert row["relative_l2_error"] < 0.002
        assert row["relative_free_residual_inf"] < 1.0e-10
        assert row["ohmic_loss"] > 0.0
        assert abs(row["center_solution"][1]) > 0.1
        assert len(row["mesh_sha256"]) == 64
        assert len(row["operator_sha256"]) == 64
        assert len(row["case_sha256"]) == 64


@pytest.mark.parametrize(
    "case, message",
    [
        ({"planar_frequency_hz": 0.0}, "positive"),
        ({"axisymmetric_frequency_hz": -1.0}, "positive"),
        ({"method": "static"}, "complex_magnetic_diffusion"),
        ({"order": 0}, "order"),
    ],
)
def test_invalid_harmonic_requests_fail_closed(case, message):
    with pytest.raises(ValueError, match=message):
        run_harmonic_magnetic_validation(case)


def test_stale_expected_result_identity_is_rejected():
    result = run_harmonic_magnetic_validation({"expected_result_sha256": "0" * 64})
    assert result["pass"] is False
    assert result["issues"] == ["expected_result_identity_matches"]
