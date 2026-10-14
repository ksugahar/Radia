import json

from radia_mcp.radia_ngsolve.field_profile_gate import symmetric_complex_field_curve_gate
from radia_mcp.radia_ngsolve.server import symmetric_complex_field_curve_gate as mcp_curve_gate


AXIS = [-3.5, -2.5, -1.5, -0.5, 0.5, 1.5, 2.5, 3.5]
FIELD_REAL = [0.15, 0.35, 0.70, 1.0, 1.0, 0.70, 0.35, 0.15]
FIELD_IMAG = [0.01, 0.02, 0.04, 0.05, 0.05, 0.04, 0.02, 0.01]


def test_symmetric_complex_curve_accepts_even_origin_bracketed_profile():
    result = symmetric_complex_field_curve_gate(
        AXIS,
        FIELD_REAL,
        FIELD_IMAG,
        log10_relative_residual=-12.3,
        min_sample_count=7,
    )
    assert result["status"] == "ok"
    assert result["checks"]["origin_sampled_or_bracketed"] is True
    assert result["metrics"]["field_symmetry_relative"] == 0.0


def test_symmetric_complex_curve_rejects_phase_asymmetry_and_weak_residual():
    bad_imag = FIELD_IMAG.copy()
    bad_imag[1] = 0.3
    result = symmetric_complex_field_curve_gate(
        AXIS,
        FIELD_REAL,
        bad_imag,
        log10_relative_residual=-4.0,
        min_sample_count=7,
    )
    assert result["status"] == "needs_attention"
    assert result["checks"]["complex_field_is_mirror_symmetric"] is False
    assert result["checks"]["solver_residual_converged"] is False


def test_symmetric_complex_curve_mcp_tool_rejects_mismatched_arrays():
    result = json.loads(mcp_curve_gate(AXIS, FIELD_REAL[:-1], -12.0))
    assert result["status"] == "invalid_input"
    assert "equal length" in result["error"]
