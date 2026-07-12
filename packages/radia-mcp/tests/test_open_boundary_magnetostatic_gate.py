from __future__ import annotations

import copy
import json

import pytest

from radia_mcp.radia_ngsolve.open_boundary_magnetostatic_gate import (
    magnetostatic_open_boundary_equivalence_gate,
)
from radia_mcp.radia_ngsolve.server import (
    magnetostatic_open_boundary_equivalence_gate as mcp_gate,
)


def _rows() -> list[dict]:
    return [
        {
            "id": "compactified_domain",
            "a": [0.0172653242, 0.0166223965, 0.0091892345, 0.0082599598, 0.0082508454],
            "bx": [-0.0395572412, -0.0243190724, 0.3003239493, -0.0159387633, -0.0149618150],
            "by": [4.22127e-5, -2.50375e-4, -4.17284e-6, 0.3649772362, -0.3654424719],
            "energy": [501.9599195, 0.2178335035],
            "coenergy": [172.9219698, 0.2177069794],
            "forces": [{"x": 0.1096563, "y": 1295.3704882}, {"x": -0.0157991, "y": -1293.9875605}],
            "mesh": {"nodes": 12264, "elements": 23924},
        },
        {
            "id": "asymptotic_shell",
            "a": [0.0219715826, 0.0213277074, 0.0138985307, 0.0129621159, 0.0129652065],
            "bx": [-0.0392945592, -0.0245084200, 0.3002515843, -0.0159426515, -0.0163265119],
            "by": [1.12847e-4, 1.12261e-4, -1.82201e-5, 0.3647576729, -0.3645446071],
            "energy": [501.9760008, 0.2178093725],
            "coenergy": [172.8828143, 0.2176829009],
            "forces": [{"x": -0.1787882, "y": 1294.7301471}, {"x": 0.2506411, "y": -1293.7704133}],
            "mesh": {"nodes": 13627, "elements": 26892},
        },
    ]


def test_gate_accepts_additive_gauge_difference_and_physical_observables() -> None:
    result = magnetostatic_open_boundary_equivalence_gate(_rows())
    assert result["status"] == "ok"
    assert result["metrics"]["max_direct_a_relative_error_diagnostic_only"] > 0.2
    assert result["checks"]["vector_potential_additive_offset_consistent"] is True
    assert result["wave_boundary_policy"] == "not_applicable_do_not_infer_from_this_gate"


def test_gate_rejects_physical_field_mismatch_despite_consistent_a_offset() -> None:
    rows = copy.deepcopy(_rows())
    rows[1]["bx"][2] *= 0.8
    result = magnetostatic_open_boundary_equivalence_gate(rows)
    assert result["status"] == "needs_attention"
    assert result["checks"]["dominant_flux_density_components_agree"] is False


def test_gate_rejects_identical_mesh_inventories() -> None:
    rows = copy.deepcopy(_rows())
    rows[1]["mesh"] = dict(rows[0]["mesh"])
    result = magnetostatic_open_boundary_equivalence_gate(rows)
    assert result["status"] == "needs_attention"
    assert result["checks"]["positive_distinct_mesh_inventories"] is False


def test_gate_refuses_wave_boundary_inference() -> None:
    with pytest.raises(ValueError, match="magnetostatic"):
        magnetostatic_open_boundary_equivalence_gate(
            _rows(), physics_regime="time_harmonic_wave"
        )


def test_mcp_tool_dispatches_the_gate() -> None:
    result = json.loads(mcp_gate(_rows()))
    assert result["status"] == "ok"
