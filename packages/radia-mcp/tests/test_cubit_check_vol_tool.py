"""Tests for the cubit_check_vol MCP tool (canonical check-vol gate).

The error path runs in any environment (structured JSON error, never a
raw traceback).  The happy path needs ngsolve + cubit-mesh-export and is
skipped in the minimal-dep matrix.
"""

import json

import pytest

from radia_mcp.cubit.server import cubit_check_vol


def test_missing_vol_returns_structured_error(tmp_path):
    out = json.loads(cubit_check_vol(str(tmp_path / "does_not_exist.vol")))
    assert out["status"] == "error"
    assert out["stage"] in {"import", "input", "check"}
    assert "error" in out


def test_explicit_missing_json_reference_is_error(tmp_path):
    pytest.importorskip("ngsolve")
    pytest.importorskip("cubit_mesh_export")
    vol = tmp_path / "cube.vol"
    _write_unit_cube_vol(vol)
    out = json.loads(cubit_check_vol(
        str(vol), json_path=str(tmp_path / "missing.json")))
    assert out["status"] == "error"
    assert out["stage"] == "input"


def _write_unit_cube_vol(vol_path):
    from netgen.occ import Box, Pnt, OCCGeometry
    from ngsolve import Mesh, TaskManager

    box = Box(Pnt(0, 0, 0), Pnt(1, 1, 1))
    box.mat("cube")
    with TaskManager():
        mesh = Mesh(OCCGeometry(box).GenerateMesh(maxh=0.5))
    mesh.ngmesh.Save(str(vol_path))


def test_unit_cube_passes_gate(tmp_path):
    pytest.importorskip("ngsolve")
    pytest.importorskip("cubit_mesh_export")
    vol = tmp_path / "cube.vol"
    _write_unit_cube_vol(vol)

    report = tmp_path / "vol_check.json"
    out = json.loads(cubit_check_vol(str(vol), report_json=str(report)))
    assert out["status"] == "ok"
    assert out["passed"] is True
    assert out["warnings"] == []
    assert out["mesh"]["n_elements"] > 0
    assert any(m["name"] == "cube" for m in out["materials"])
    # report artifact written and self-consistent
    assert report.is_file()
    payload = json.loads(report.read_text(encoding="utf-8"))
    assert payload["passed"] is True
