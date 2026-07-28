from __future__ import annotations

import pytest

from radia_mcp.radia_ngsolve.age_retirement_validation import run_age_retirement_validation


@pytest.fixture(scope="module")
def verified():
    return run_age_retirement_validation({})


def test_generated_vol_executes_no_remesh_age_sweep(verified):
    assert verified["status"] == "verified"
    assert verified["validated_capabilities"] == ["moving_band", "periodic_boundary"]
    assert verified["torque_summary"]["closure_relative_error"] < 1.0e-8
    assert verified["torque_summary"]["phase_sign_reversal_observed"] is True


def test_execution_identities_and_reuse_are_explicit(verified):
    for key in ("operator_sha256", "age_factorization_sha256", "angle_grid_sha256", "excitation_sha256", "torque_output_sha256", "result_sha256"):
        assert len(verified[key]) == 64
    assert verified["checks"]["mesh_operator_and_factorization_reused"] is True


@pytest.mark.parametrize(
    "case, message",
    [
        ({"method": "remesh_each_angle"}, "ngsolve_age_phase_only"),
        ({"angle_samples": 3}, "angle_samples"),
        ({"boundary": "periodic"}, "periodic boundary sign"),
        ({"boundary_phase": 1}, "boundary_phase"),
    ],
)
def test_invalid_motion_contract_fails_closed(case, message):
    with pytest.raises(ValueError, match=message):
        run_age_retirement_validation(case)


def test_stale_result_identity_is_rejected():
    result = run_age_retirement_validation({"expected_result_sha256": "0" * 64})
    assert result["pass"] is False
    assert result["issues"] == ["expected_result_identity_matches"]
