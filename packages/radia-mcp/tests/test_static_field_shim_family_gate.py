from __future__ import annotations

import copy
import json

from radia_mcp.radia_ngsolve.server import static_field_shim_family_gate as mcp_gate
from radia_mcp.radia_ngsolve.static_field_shim_family_gate import (
    static_field_shim_family_gate,
)


def _case(role: str, center_z: float, p2p_ppm: float, divergence: float) -> dict:
    return {
        "role": role,
        "center_b_t": [1.0e-7, -2.0e-7, center_z],
        "central_axial_peak_to_peak_ppm": p2p_ppm,
        "central_axial_rms_ppm": 0.3 * p2p_ppm,
        "center_transverse_relative": 3.0e-6,
        "central_divergence_max_relative": divergence,
        "grid_shape": [21, 21, 20] if role == "paired_source" else [21, 21, 21],
        "row_count": 8820 if role == "paired_source" else 9261,
        "central_sample_count": 27,
        "roi_half_width_mm": 10.0,
        "coordinate_unit": "mm",
        "field_unit": "T",
    }


def _cases() -> list[dict]:
    return [
        _case("single_source", 0.10, 100000.0, 0.01),
        _case("paired_source", 0.16, 30000.0, 0.008),
        _case("balanced_shim", 0.09, 40000.0, 0.02),
        _case("offset_shim", 0.085, 35000.0, 0.021),
    ]


def test_static_field_shim_family_gate_accepts_ordering_and_map_quality() -> None:
    result = static_field_shim_family_gate(_cases())
    assert result["status"] == "ok"
    assert result["metrics"]["paired_source_field_ratio"] > 1.5


def test_static_field_shim_family_gate_rejects_lost_uniformity_gain() -> None:
    cases = copy.deepcopy(_cases())
    cases[1]["central_axial_peak_to_peak_ppm"] = 60000.0
    cases[1]["central_axial_rms_ppm"] = 20000.0
    result = static_field_shim_family_gate(cases)
    assert result["status"] == "needs_attention"
    assert result["checks"]["paired_source_improves_roi_uniformity"] is False


def test_static_field_shim_family_gate_rejects_unresolved_shim_pair() -> None:
    cases = copy.deepcopy(_cases())
    cases[3]["center_b_t"][2] = cases[2]["center_b_t"][2]
    result = static_field_shim_family_gate(cases)
    assert result["status"] == "needs_attention"
    assert result["checks"]["shim_variants_are_resolved"] is False


def test_static_field_shim_family_mcp_tool_dispatches() -> None:
    result = json.loads(mcp_gate(_cases()))
    assert result["status"] == "ok"
