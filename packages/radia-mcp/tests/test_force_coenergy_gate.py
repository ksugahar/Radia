import json
import math

import pytest

from radia_mcp.radia_ngsolve.force_coenergy_gate import force_coenergy_displacement_gate
from radia_mcp.radia_ngsolve.server import force_coenergy_displacement_gate as mcp_gate


def _quadratic_case():
    positions = [0.002 * index for index in range(7)]
    coenergy = [2.0 - 40.0 * x + 100.0 * x * x for x in positions]
    forces = [-40.0 + 200.0 * x for x in positions]
    return positions, coenergy, forces


def test_force_coenergy_gate_accepts_constant_current_virtual_work_identity():
    positions, coenergy, forces = _quadratic_case()
    result = force_coenergy_displacement_gate(positions, coenergy, forces)
    assert result["status"] == "ok"
    assert result["max_central_relative_error"] < 1.0e-12
    assert result["endpoint_errors_are_diagnostic_only"] is True
    assert result["rows"][0]["stencil"] == "forward"


def test_force_coenergy_gate_rejects_force_with_wrong_projection_sign():
    positions, coenergy, forces = _quadratic_case()
    result = force_coenergy_displacement_gate(positions, coenergy, [-f for f in forces])
    assert result["status"] == "needs_attention"
    assert result["checks"]["central_virtual_work_matches_direct_force"] is False


def test_force_coenergy_mcp_tool_dispatches_json_and_handles_bad_shape():
    positions, coenergy, forces = _quadratic_case()
    result = json.loads(mcp_gate(positions, coenergy, forces))
    assert result["status"] == "ok"
    bad = json.loads(mcp_gate(positions, coenergy[:-1], forces))
    assert bad["status"] == "invalid_input"


def test_force_coenergy_gate_requires_constant_current_semantics():
    positions, coenergy, forces = _quadratic_case()
    result = force_coenergy_displacement_gate(
        positions, coenergy, forces, energy_kind="stored_energy_at_fixed_flux"
    )
    assert result["status"] == "needs_attention"
    assert result["checks"]["constant_current_coenergy_recorded"] is False


@pytest.mark.parametrize(
    "case_id",
    [
        "v7_public_force_energy_derivative_sign_conflict",
        "v7_public_axisymmetric_two_pi_double_count",
    ],
)
def test_generalization_v7_public(case_id):
    positions, coenergy, forces = _quadratic_case()
    if case_id == "v7_public_force_energy_derivative_sign_conflict":
        forces = [-force for force in forces]
    else:
        forces = [2.0 * math.pi * force for force in forces]
    result = force_coenergy_displacement_gate(positions, coenergy, forces)
    assert result["status"] == "needs_attention"
    assert result["checks"]["central_virtual_work_matches_direct_force"] is False
