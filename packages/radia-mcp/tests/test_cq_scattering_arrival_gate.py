import json

from radia_mcp.radia_ngsolve.cq_scattering_arrival_gate import cq_scattering_arrival_gate
from radia_mcp.radia_ngsolve.server import cq_scattering_arrival_gate as mcp_gate


def _args():
    return dict(time_step_s=0.32, geometric_arrival_s=4.0, measured_peak_s=4.48, max_relative_residual=2e-16, finite_response=True, real_time_response=True)


def test_cq_arrival_gate_accepts_causal_peak_near_ray_time():
    result = cq_scattering_arrival_gate(**_args())
    assert result["status"] == "ok"
    assert json.loads(mcp_gate(**_args()))["status"] == "ok"


def test_cq_arrival_gate_rejects_acausal_peak():
    result = cq_scattering_arrival_gate(**{**_args(), "measured_peak_s": 3.2})
    assert result["status"] == "needs_attention"
    assert result["checks"]["peak_not_acausal"] is False


def test_cq_arrival_gate_rejects_complex_or_unconverged_response():
    result = cq_scattering_arrival_gate(**{**_args(), "real_time_response": False, "max_relative_residual": 1e-3})
    assert result["checks"]["inverse_transform_real"] is False
    assert result["checks"]["linear_solves_converged"] is False
