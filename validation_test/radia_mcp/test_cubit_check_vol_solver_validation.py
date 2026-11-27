"""Solver-backed validation for the canonical cubit_check_vol gate."""

from __future__ import annotations

import json

import pytest

from radia_mcp.cubit.server import cubit_check_vol


def _write_unit_cube_vol(vol_path):
    from netgen.occ import Box, OCCGeometry, Pnt
    from ngsolve import Mesh, TaskManager

    box = Box(Pnt(0, 0, 0), Pnt(1, 1, 1))
    box.mat("cube")
    with TaskManager():
        mesh = Mesh(OCCGeometry(box).GenerateMesh(maxh=0.5))
    mesh.ngmesh.Save(str(vol_path))


def test_explicit_missing_json_reference_is_error(tmp_path):
    pytest.importorskip("ngsolve")
    pytest.importorskip("cubit_mesh_export")
    vol = tmp_path / "cube.vol"
    _write_unit_cube_vol(vol)
    out = json.loads(
        cubit_check_vol(str(vol), json_path=str(tmp_path / "missing.json"))
    )
    assert out["status"] == "error"
    assert out["stage"] == "input"


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
    assert any(material["name"] == "cube" for material in out["materials"])
    assert report.is_file()
    payload = json.loads(report.read_text(encoding="utf-8"))
    assert payload["passed"] is True
