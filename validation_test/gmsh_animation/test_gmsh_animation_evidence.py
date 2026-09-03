from __future__ import annotations

import json
import math
from pathlib import Path


HERE = Path(__file__).resolve().parent


def _load(name: str) -> dict:
    return json.loads((HERE / name).read_text(encoding="utf-8"))


def test_inspection_evidence_preserves_mesh_and_displacement_contract():
    result = _load("gmsh_animation_results.json")
    assert result["schema"] == "radia.validation.gmsh_animation.inspect.v1"
    assert result["msh"]["mesh_format"] == "4.1 0 8"
    assert result["msh"]["node_data_count"] == 21
    assert all(result["checks"].values())
    assert math.isclose(result["checks"]["final_displacement_m"], 0.15)


def test_export_evidence_preserves_sidecar_and_media_contract():
    result = _load("gmsh_animation_export_results.json")
    assert result["schema"] == "radia.validation.gmsh_animation_export.v1"
    assert result["animation"]["node_data_frames"] == 21
    assert result["animation"]["exported_frame_count"] == 21
    assert all(result["checks"].values())
    assert all(row["exists"] and row["bytes"] > 0 for row in result["artifacts"].values())
    assert all(value >= 0.0 for value in result["timing_breakdown_s"].values())
