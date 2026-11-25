"""Guards: native AGE rows round-trip through both existing comparison contracts.

Assumes: uniformly sampled, unskewed, zero-current synthetic traces.
"""

import json
from math import cos, sin, tau

import pytest

from radia_mcp.motor.back_emf_cross_lane_gate import back_emf_cross_lane_gate
from radia_mcp.motor.cogging_torque_cross_lane_gate import cogging_torque_cross_lane_gate
from radia_mcp.motor.server import motor_age_lanes

EMF_LIMITS = dict(grid_deg_e=1e-6, amplitude_relative=0.01, phase_deg_e=1,
                  rms_relative=0.01, thd_absolute=0.01)
TORQUE_LIMITS = dict(peak_to_peak_nm=1e-6, waveform_rms_nm=1e-6, harmonic_nm=1e-6)


def test_age_export_uses_independent_periods_and_preserves_sign():
    """Guards: both exported contracts and their physical metadata.
    Assumes: a four-pole six-slot model with declared torque symmetry.
    """
    n = 24
    back = [dict(theta_r=tau*i/(2*n), converged=True,
                 linkage_stack_basis="total", phase_flux_wb=[0.02*cos(tau*(i/n-j/3)) for j in range(3)]) for i in range(n)]
    torque = [dict(theta_r=tau*i/(12*n), converged=True, axial_length_m=0.05, torque=-0.02*sin(tau*i/n))
              for i in range(n)]
    payload = dict(back_emf_rows=back, cogging_rows=torque, stack_length_m=0.05, pole_pairs=2,
                   phase_current_A=0, speed_rpm=1200, stack_basis="total",
                   linkage_polarity="positive_winding", slots=6, slot_pole_symmetry=True,
                   torque_basis="total", model_fraction=1)
    result = json.loads(motor_age_lanes(json.dumps(payload)))
    linkage, cogging = result["back_emf_lane"], result["cogging_lane"]
    assert linkage["source_angle_unit"] == "rad"
    assert linkage["mechanical_angles_deg"][-1] == pytest.approx(180*23/24)
    assert cogging["mechanical_angles_deg"][-1] == pytest.approx(30*23/24)
    assert cogging["torque_nm"][6] == pytest.approx(-0.02)
    assert back_emf_cross_lane_gate(linkage, linkage, EMF_LIMITS)["status"] == "pass"
    assert cogging_torque_cross_lane_gate(cogging, cogging, TORQUE_LIMITS)["status"] == "pass"
    inverted = {**linkage, "phase_flux_wb": [[-x for x in row] for row in linkage["phase_flux_wb"]]}
    assert back_emf_cross_lane_gate(linkage, inverted, EMF_LIMITS)["status"] == "mismatch"


def test_slotless_samples_stay_slotless_and_invalid_sweeps_are_rejected():
    """Guards: slotless traces retain their physical identity and failed solves do not export.
    Assumes: the cogging gate requires slots greater than zero.
    """
    payload = dict(cogging_rows=[dict(theta_r=0, torque=-0.001, converged=True, axial_length_m=1.0)],
                   stack_length_m=0.05, pole_pairs=2, phase_current_A=0, slots=0,
                   slot_pole_symmetry=False, torque_basis="per_metre", model_fraction=1)
    result = json.loads(motor_age_lanes(json.dumps(payload)))
    assert result["cogging_unavailable"]["torque_nm"] == [-0.001]
    payload["phase_current_A"] = 1
    assert json.loads(motor_age_lanes(json.dumps(payload)))["status"] == "invalid_input"
    payload["phase_current_A"] = 0
    payload["cogging_rows"][0]["converged"] = False
    assert json.loads(motor_age_lanes(json.dumps(payload)))["status"] == "invalid_input"


@pytest.mark.parametrize(("rows", "error"), [
    ([None], "row 0 must be an object"),
    ([{"theta_r": 0, "converged": True}], "row 0 is missing torque"),
    ([{"theta_r": 0, "converged": True, "torque": 0, "axial_length_m": 0.1}],
     "row 0 axial_length_m does not match torque_basis"),
    ([{"theta_r": 0, "converged": True, "torque": 0, "axial_length_m": 1.0}] * 2,
     "row 1 angle must be strictly increasing"),
])
def test_invalid_cogging_rows_report_the_failure(rows, error):
    payload = dict(cogging_rows=rows, stack_length_m=0.05, pole_pairs=2,
                   phase_current_A=0, slots=6, slot_pole_symmetry=True,
                   torque_basis="per_metre", model_fraction=1)
    result = json.loads(motor_age_lanes(json.dumps(payload)))
    assert result["status"] == "invalid_input"
    assert error in result["error"]


def test_linkage_basis_must_match():
    payload = dict(back_emf_rows=[dict(theta_r=0, converged=True, phase_flux_wb=[1, 2, 3],
                                      linkage_stack_basis="per_metre")],
                   stack_length_m=0.05, pole_pairs=2, phase_current_A=0, speed_rpm=1200,
                   stack_basis="total", linkage_polarity="positive_winding")
    result = json.loads(motor_age_lanes(json.dumps(payload)))
    assert result["status"] == "invalid_input"
    assert "row 0 linkage_stack_basis does not match" in result["error"]


def test_age_sweep_stub_preserves_default_rows_and_exports_provenance(monkeypatch):
    from radia_mcp.radia_ngsolve import airgap_motor_workflow as workflow

    sources = []
    field = object()

    def solve(*args, source_lf, **kwargs):
        sources.append(source_lf)
        return field, {"niter_used": 2, "converged": True}

    monkeypatch.setattr(workflow, "age_motor_nonlinear_solve", solve)
    monkeypatch.setattr(workflow, "airgap_torque", lambda coupling, solution, axial_length: -axial_length)
    args = (None, {}, None, lambda theta: None, [0.1])
    old = workflow.age_motor_rotation_sweep(*args)
    assert old == [dict(theta_r=0.1, torque=-1, niter_used=2, converged=True)]
    exported = workflow.age_motor_rotation_sweep(
        *args, axial_length=0.05, source_fn=lambda theta: ("source", theta),
        phase_linkage_fn=lambda solution, theta: (theta, -theta, 0),
        linkage_stack_basis="total", record_export_metadata=True,
    )
    assert sources == [None, ("source", 0.1)]
    assert exported == [dict(theta_r=0.1, torque=-0.05, niter_used=2, converged=True,
                             axial_length_m=0.05, phase_flux_wb=(0.1, -0.1, 0),
                             linkage_stack_basis="total")]


def test_age_sweep_aborts_when_stub_solve_raises(monkeypatch):
    from radia_mcp.radia_ngsolve import airgap_motor_workflow as workflow

    angles = []

    def fail(*args, source_lf, **kwargs):
        angles.append(source_lf)
        raise RuntimeError("Picard did not converge")

    monkeypatch.setattr(workflow, "age_motor_nonlinear_solve", fail)
    with pytest.raises(RuntimeError, match="Picard did not converge"):
        workflow.age_motor_rotation_sweep(
            None, {}, None, lambda theta: None, [0.1, 0.2], source_fn=lambda theta: theta,
        )
    assert angles == [0.1]
