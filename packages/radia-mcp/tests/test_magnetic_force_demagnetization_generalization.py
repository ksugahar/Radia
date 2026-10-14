"""Magnetic force / demagnetization artifact-identity gate (v45, v46 roots).

One positive closure per independent payload root; every negative test keeps
its own failure signal.  Payload builders live in
``_magnetic_force_demagnetization_generalization_payloads``.
"""

from __future__ import annotations

from _magnetic_force_demagnetization_generalization_payloads import _identity_v45, _identity_v46
from radia_mcp.radia_ngsolve.magnetic_force_v44_identity import validate_public_identity


def test_v45_magnetic_force_public_identity_accepts_closed_artifacts():
    checks = validate_public_identity(_identity_v45())
    assert checks and all(checks.values())


def test_v45_magnetic_force_public_identity_rejects_dynamic_phase_mutation():
    identity = _identity_v45()
    identity["v45_public_magnetic_bearing_dynamic_stiffness_phase_damping_force_power_stability_owner_mismatch"]["result_phase_deg"] = [0.0, -10.0]
    checks = validate_public_identity(identity)
    assert checks and not all(checks.values())


def test_v46_public_identity_accepts_closed_artifacts():
    checks = validate_public_identity(_identity_v46())
    assert checks and all(checks.values())


def test_v46_public_identity_rejects_partial_nan_and_restart_mutations():
    identity = _identity_v46()
    identity["v46_public_magnetic_force_partial_solve_unit_scale_coordinate_frame_nan_mismatch"]["result_finite_values"] = False
    identity["v46_public_demagnetization_curve_branch_restart_temperature_window_mismatch"]["result_partial_path_status"] = "partial"
    checks = validate_public_identity(identity)
    assert checks and not all(checks.values())
