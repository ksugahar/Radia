import copy
import json

from radia_mcp.build123d.platonic_family_gate import (
    EXPECTED,
    platonic_solid_family_gate,
)
from radia_mcp.build123d.server import build123d_platonic_solid_family_gate as mcp_gate


def _summary():
    rows = []
    for name, (vertices, edges, faces) in EXPECTED.items():
        rows.append(
            {
                "name": name,
                "vertices": vertices,
                "edges": edges,
                "faces": faces,
                "analytic_relative_error": 2.0e-16,
                "self_roundtrip_relative_error": 2.0e-12,
                "external_relative_error": 2.0e-12,
                "external_volume_count": 1,
                "valid": True,
            }
        )
    return {
        "upstream_commit": "a" * 40,
        "source_sha256": "b" * 64,
        "source_execution_mode": "exact_source_with_display_stub",
        "shape_valid_access": "property",
        "rows": rows,
    }


def test_accepts_complete_platonic_family_contract():
    result = platonic_solid_family_gate(_summary())
    assert result["status"] == "ok"
    assert json.loads(mcp_gate(json.dumps(_summary())))["status"] == "ok"


def test_rejects_equal_volume_with_wrong_topology():
    bad = copy.deepcopy(_summary())
    bad["rows"][-1]["faces"] = 19
    result = platonic_solid_family_gate(bad)
    assert result["status"] == "needs_attention"
    assert set(result["issues"]) >= {
        "known_vertex_edge_face_counts",
        "euler_characteristic_two",
    }


def test_rejects_external_translation_loss_and_wrong_api_contract():
    bad = copy.deepcopy(_summary())
    bad["rows"][0]["external_relative_error"] = 0.1
    bad["shape_valid_access"] = "method"
    result = platonic_solid_family_gate(bad)
    assert set(result["issues"]) >= {
        "external_kernel_volumes_match",
        "installed_api_contract_recorded",
    }
