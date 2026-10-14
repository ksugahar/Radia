import json

from radia_mcp.radia_ngsolve.parallel_wire_force_refinement_gate import (
    parallel_wire_force_refinement_gate,
)
from radia_mcp.radia_ngsolve.server import (
    parallel_wire_force_refinement_gate as mcp_gate,
)


LEVELS = [1, 2, 3, 4, 5]
F1 = [
    [-1.9500468e-7, -8.7304e-10],
    [-1.9676452e-7, -1.1667e-9],
    [-1.9930409e-7, -1.3409e-10],
    [-1.9846973e-7, 6.8247e-11],
    [-1.9976838e-7, -8.5147e-11],
]
F2 = [
    [1.9542735e-7, -6.6129e-10],
    [1.9666364e-7, 6.9491e-10],
    [1.9884370e-7, 3.9145e-11],
    [1.9832165e-7, 5.1803e-11],
    [1.9959257e-7, -9.4344e-11],
]


def test_parallel_wire_force_refinement_accepts_nonmonotone_convergence():
    result = parallel_wire_force_refinement_gate(
        LEVELS,
        F1,
        F2,
        expected_force_magnitude=2.0e-7,
        expected_wire2_radial_sign=1,
    )
    assert result["status"] == "ok"
    assert result["error_is_monotone"] is False
    assert result["final_relative_error"] < 0.003
    assert result["initial_to_final_error_ratio"] > 10.0


def test_parallel_wire_force_refinement_rejects_stale_single_body_result():
    bad_pair = [row[:] for row in F1]
    result = parallel_wire_force_refinement_gate(
        LEVELS,
        F1,
        bad_pair,
        expected_force_magnitude=2.0e-7,
        expected_wire2_radial_sign=1,
    )
    assert result["status"] == "needs_attention"
    assert result["checks"]["opposite_pair_force_sign"] is False
    assert result["checks"]["final_action_reaction_residual_ok"] is False


def test_parallel_wire_force_refinement_mcp_dispatches_json():
    result = json.loads(mcp_gate(LEVELS, F1, F2, 2.0e-7, [1.0, 0.0], 1))
    assert result["status"] == "ok"
    assert result["policy"] == "parallel_wire_force_refinement_gate_v1"
