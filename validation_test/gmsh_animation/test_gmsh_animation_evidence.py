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
    # This is the archived calculation, not today's format acceptance oracle.
    # Keep its original provenance; the executable inspector is checked below.
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


def test_current_inspector_reads_blocked_displacement_data(tmp_path):
    import importlib.util

    path = HERE.parents[1] / 'docs/gmsh_animation/gmsh_animation_inspect.py'
    spec = importlib.util.spec_from_file_location('animation_inspection', path)
    inspector = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(inspector)
    msh = tmp_path / 'displacement.msh'
    msh.write_text(
        '$MeshFormat\n4.2 0 8\n$EndMeshFormat\n'
        '$NodeData\n1\n"Displacement [m]"\n1\n1\n3\n1\n3\n2\n'
        '1 2 17 19\n2\n17\n19\n0 0 0\n0.15 0 0\n$EndNodeData\n',
        encoding='utf-8')
    result = inspector.inspect_msh(msh)
    assert result['mesh_format'] == '4.2 0 8'
    assert result['node_data'][0]['max_norm'] == 0.15
    assert result['node_data'][0]['mean_norm'] == 0.075
