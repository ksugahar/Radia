from __future__ import annotations

import copy

import pytest

from radia_mcp.radia_ngsolve.vol2d_postprocess import (
    POSTPROCESS_SCHEMA,
    _canonical,
    _coordinate_contract,
    _named_paths,
    _named_points,
    _normalize_current_rows,
    _sha,
    analyze_vol2d_postprocess,
    postprocess_replay_gate,
)


def _entry(content: str, filename: str) -> dict:
    return {"content": content, "filename": filename, "sha256": _sha(content)}


def _artifact() -> dict:
    request_contract = {
        "schema": "radia.vol2d-postprocess-request.v1",
        "export_basename": "field_post",
    }
    states = [[0.25, -0.5]]
    rows = [{"row_index": 0, "branch_current_a": [1.0]}]
    contents = {
        "csv": "row_index,current_a:winding\n0,1\n",
        "gmsh_msh": "$MeshFormat\n4.1 0 8\n$EndMeshFormat\n",
        "gmsh_geo": 'Merge "field_post.msh";\n',
        "gmsh_geo_opt": "General.InitialModule = 5;\n",
        "gmsh_msh_opt": "General.InitialModule = 1;\n",
    }
    contract = {
        "schema": POSTPROCESS_SCHEMA,
        "status": "solved",
        "operation": "solve",
        "request_contract": request_contract,
        "request_contract_sha256": _sha(request_contract),
        "field_state_rows": states,
        "field_state_table_sha256": _sha(states),
        "sweep_rows": rows,
        "result_table_sha256": _sha(rows),
        "factorization_count": 1,
        "mesh_rebuild_count": 0,
        "operator_rebuild_count": 0,
        "solve_count": 1,
        "export_content_sha256": {name: _sha(content) for name, content in contents.items()},
    }
    json_content = _canonical(contract)
    exports = {
        "json": _entry(json_content, "field_post.json"),
        "csv": _entry(contents["csv"], "field_post.csv"),
        "gmsh_msh": _entry(contents["gmsh_msh"], "field_post.msh"),
        "gmsh_geo": _entry(contents["gmsh_geo"], "field_post.geo"),
        "gmsh_geo_opt": _entry(contents["gmsh_geo_opt"], "field_post.geo.opt"),
        "gmsh_msh_opt": _entry(contents["gmsh_msh_opt"], "field_post.msh.opt"),
    }
    contract["canonical_json_sha256"] = exports["json"]["sha256"]
    return {
        "schema": POSTPROCESS_SCHEMA,
        "status": "solved",
        "operation": "solve",
        "result_contract": contract,
        "exports": exports,
    }


def test_replay_accepts_complete_deterministic_artifact() -> None:
    gate = postprocess_replay_gate(_artifact())
    assert gate["status"] == "accepted"
    assert gate["pass"] is True
    assert all(gate["checks"].values())


def test_replay_rejects_tampered_export_and_state() -> None:
    artifact = copy.deepcopy(_artifact())
    artifact["result_contract"]["field_state_rows"][0][0] += 1.0
    artifact["exports"]["csv"]["content"] += "tampered\n"
    gate = postprocess_replay_gate(artifact)
    assert gate["status"] == "rejected"
    assert gate["checks"]["field_state_table_sha256"] is False
    assert gate["checks"]["csv_content_sha256"] is False


def test_current_rows_reject_duplicate_and_nonfinite_operating_points() -> None:
    with pytest.raises(ValueError, match="duplicate"):
        _normalize_current_rows([[1.0], [1.0]], 1)
    with pytest.raises(ValueError, match="finite"):
        _normalize_current_rows([[float("nan")]], 1)


def test_path_contract_preserves_order_and_rejects_degenerate_segments() -> None:
    path = _named_paths(
        [{"name": "p", "points_m": [[0.0, 0.0], [1.0, 0.0]], "samples_per_segment": 3}],
        ("x", "y"),
    )[0]
    assert path["points_m"] == [[0.0, 0.0], [1.0, 0.0]]
    assert path["closed"] is False
    with pytest.raises(ValueError, match="duplicate"):
        _named_paths(
            [{"name": "p", "points_m": [[0.0, 0.0], [0.0, 0.0]]}],
            ("x", "y"),
        )


def test_point_and_coordinate_contracts_are_explicit() -> None:
    points = _named_points(
        [{"name": "probe", "coordinates_m": [0.2, 0.3]}], ("r", "z")
    )
    assert points[0]["coordinate_names"] == ["r", "z"]
    assert _coordinate_contract("planar")["field_components"] == ["B_x", "B_y"]
    assert _coordinate_contract("axisymmetric_henrotte")["field_components"] == [
        "B_r",
        "B_z",
    ]


def test_closed_world_dispatch_rejects_unknown_operation_without_solver() -> None:
    with pytest.raises(ValueError, match="solve or replay_gate"):
        analyze_vol2d_postprocess({"operation": "arbitrary"})
