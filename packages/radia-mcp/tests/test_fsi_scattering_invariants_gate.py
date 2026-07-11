import json

from radia_mcp.radia_ngsolve.fsi_scattering_invariants_gate import (
    fsi_scattering_invariants_gate,
)
from radia_mcp.radia_ngsolve.server import fsi_scattering_invariants_gate as mcp_gate


def _args():
    return {
        "reciprocity_relative_error": 6.4e-4,
        "optical_theorem_relative_error": 8.2e-4,
        "bem_dtn_relative_error": 1.32e-2,
        "max_solver_residual": 8.6e-15,
        "lossless_material": True,
        "time_convention": "exp(-i omega t)",
    }


def test_gate_accepts_lossless_reciprocal_energy_closure():
    result = fsi_scattering_invariants_gate(**_args())
    assert result["status"] == "ok"
    assert json.loads(mcp_gate(**_args()))["status"] == "ok"


def test_gate_rejects_optical_theorem_or_exterior_disagreement():
    result = fsi_scattering_invariants_gate(
        **{**_args(), "optical_theorem_relative_error": 0.08, "bem_dtn_relative_error": 0.07}
    )
    assert result["status"] == "needs_attention"
    assert result["checks"]["optical_theorem_energy_closure"] is False
    assert result["checks"]["p1_bem_high_order_dtn_agreement"] is False


def test_gate_rejects_absorbing_or_ambiguous_convention_claim():
    result = fsi_scattering_invariants_gate(
        **{**_args(), "lossless_material": False, "time_convention": "unspecified"}
    )
    assert result["checks"]["lossless_material_declared"] is False
    assert result["checks"]["outgoing_time_convention_explicit"] is False
