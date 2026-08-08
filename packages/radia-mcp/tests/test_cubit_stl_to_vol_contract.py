"""cubit_stl_to_vol input contract + the gmsh volume referee.

No Cubit needed: every test here exercises the fail-fast input paths and
the gmsh-side total-volume integrator the closure gate relies on.  The
live meshing path runs in validation_test (needs a Cubit license).
"""
from __future__ import annotations

import json
from pathlib import Path

import pytest

pytest.importorskip("gmsh")

from radia_mcp.cubit import server  # noqa: E402
from radia_mcp.cubit.stl_inspect import inspect_stl  # noqa: E402
from radia_mcp.gmsh.msh_inspect import mesh_total_volume  # noqa: E402


def _call(fn, **kw):
    r = fn(**kw)
    return json.loads(r) if isinstance(r, str) else r


def _write_ascii_stl(path: Path, triangles) -> None:
    rows = ["solid fixture"]
    for triangle in triangles:
        rows.extend(["  facet normal 0 0 0", "    outer loop"])
        rows.extend(f"      vertex {x} {y} {z}" for x, y, z in triangle)
        rows.extend(["    endloop", "  endfacet"])
    rows.append("endsolid fixture")
    path.write_text("\n".join(rows) + "\n", encoding="ascii")


def _write_unit_cube_stl(path: Path) -> None:
    v = [(0, 0, 0), (1, 0, 0), (1, 1, 0), (0, 1, 0),
         (0, 0, 1), (1, 0, 1), (1, 1, 1), (0, 1, 1)]
    faces = [(0, 2, 1), (0, 3, 2), (4, 5, 6), (4, 6, 7),
             (0, 1, 5), (0, 5, 4), (3, 7, 6), (3, 6, 2),
             (0, 4, 7), (0, 7, 3), (1, 2, 6), (1, 6, 5)]
    _write_ascii_stl(path, [[v[index] for index in face] for face in faces])


def test_missing_stl_is_input_error(tmp_path):
    r = _call(server.cubit_stl_to_vol, stl_path=str(tmp_path / "nope.stl"))
    assert r["status"] == "error" and r["kind"] == "input"


def test_bad_scheme_is_input_error(tmp_path):
    p = tmp_path / "x.stl"
    p.write_text("solid x\nendsolid x\n")
    r = _call(server.cubit_stl_to_vol, stl_path=str(p), scheme="prism")
    assert r["status"] == "error" and r["kind"] == "input"
    assert "hex" in r["error"] and "tet" in r["error"]


def test_open_surface_is_input_error(tmp_path):
    p = tmp_path / "open.stl"
    _write_ascii_stl(p, [[(0, 0, 0), (1, 0, 0), (0, 1, 0)]])
    r = _call(server.cubit_stl_to_vol, stl_path=str(p))
    assert r["status"] == "error" and r["kind"] == "input"
    assert "watertight" in r["error"]


def test_stl_inspector_recovers_closed_unit_cube_volume(tmp_path):
    p = tmp_path / "cube.stl"
    _write_unit_cube_stl(p)

    r = inspect_stl(p)

    assert r["watertight"] is True
    assert r["winding_consistent"] is True
    assert r["triangle_count"] == 12
    assert r["open_edge_count"] == 0
    assert r["nonmanifold_edge_count"] == 0
    assert abs(r["volume"] - 1.0) < 1e-12


def test_stl_route_uses_plugin_aware_headless_process(tmp_path, monkeypatch):
    p = tmp_path / "cube.stl"
    _write_unit_cube_stl(p)
    seen = {}

    def fake_headless(commands, **kwargs):
        seen["commands"] = list(commands)
        seen["kwargs"] = kwargs
        return {"status": "completed", "exit_code": 0}

    monkeypatch.setattr(server._cs, "run_headless_journal", fake_headless)
    r = _call(server.cubit_stl_to_vol, stl_path=str(p),
              out_vol=str(tmp_path / "cube.vol"),
              out_msh=str(tmp_path / "cube.msh"))

    assert r["status"] == "error" and r["stage"] == "cubit"
    assert any(line.startswith("export netgen") for line in seen["commands"])
    assert any("overwrite" in line for line in seen["commands"])
    assert seen["kwargs"]["working_directory"] == tmp_path


def _write_unit_cube_msh(path: Path) -> None:
    import gmsh

    gmsh.initialize(["-noconfig"])
    try:
        gmsh.option.setNumber("General.Terminal", 0)
        gmsh.model.add("cube")
        gmsh.model.occ.addBox(0, 0, 0, 1, 1, 1)
        gmsh.model.occ.synchronize()
        gmsh.option.setNumber("Mesh.MeshSizeMin", 0.5)
        gmsh.option.setNumber("Mesh.MeshSizeMax", 0.5)
        gmsh.model.mesh.generate(3)
        gmsh.write(str(path))
    finally:
        gmsh.finalize()


def test_mesh_total_volume_unit_cube(tmp_path):
    p = tmp_path / "cube.msh"
    _write_unit_cube_msh(p)
    r = mesh_total_volume(p)
    assert r["ran"] and r["ok"], r
    assert abs(r["total_volume"] - 1.0) < 1e-9
    assert r["n_elements_3d"] > 0
    assert r["min_jacobian_det"] > 0.0


def test_mesh_total_volume_missing_file(tmp_path):
    r = mesh_total_volume(tmp_path / "none.msh")
    assert r["ok"] is False and "not found" in r["error"]


def test_vol_surface_element_count_reads_section(tmp_path):
    p = tmp_path / "toy.vol"
    p.write_text(
        "mesh3d\n\nvolumeelements\n1\n1 4 1 2 3 4\n\n"
        "surfaceelements\n3\nrow\nrow\nrow\n",
        encoding="utf-8")
    assert server._vol_surface_element_count(p) == 3


def test_vol_surface_element_count_missing_section_is_negative(tmp_path):
    p = tmp_path / "bald.vol"
    p.write_text("mesh3d\n\nvolumeelements\n1\n1 4 1 2 3 4\n",
                 encoding="utf-8")
    assert server._vol_surface_element_count(p) == -1


def test_hex_route_preserves_free_sideset_and_binds_mesh_geometry(
        tmp_path, monkeypatch):
    """Keep the Sculpt skin as both a free sideset and mesh geometry.

    Older exporters ignored free faces and produced a surface-less `.vol`,
    making a charge-based solve silently demag-free.
    """
    p = tmp_path / "cube.stl"
    _write_unit_cube_stl(p)
    seen = {}

    def fake_headless(commands, **kwargs):
        seen["commands"] = list(commands)
        return {"status": "completed", "exit_code": 0}

    monkeypatch.setattr(server._cs, "run_headless_journal", fake_headless)
    r = _call(server.cubit_stl_to_vol, stl_path=str(p), scheme="hex",
              out_vol=str(tmp_path / "cube.vol"),
              out_msh=str(tmp_path / "cube.msh"))

    assert r["status"] == "error" and r["stage"] == "cubit"
    commands = seen["commands"]
    sculpt = next(i for i, line in enumerate(commands)
                  if line.startswith("sculpt volume all"))
    assert "gen_sidesets 2" in commands[sculpt]
    assert commands[sculpt + 1] == \
        "create mesh geometry hex all feature_angle 135"
    assert commands[sculpt + 2] == "delete volume with not is_meshed"
