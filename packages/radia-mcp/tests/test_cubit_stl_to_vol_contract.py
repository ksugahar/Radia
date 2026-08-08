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

from radia_mcp.cubit import server
from radia_mcp.cubit.stl_inspect import inspect_stl
from radia_mcp.gmsh.msh_inspect import mesh_total_volume


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


def _unit_cube_triangles(*, origin=(0.0, 0.0, 0.0), size=1.0,
                         reverse=False):
    v = [(0, 0, 0), (1, 0, 0), (1, 1, 0), (0, 1, 0),
         (0, 0, 1), (1, 0, 1), (1, 1, 1), (0, 1, 1)]
    v = [(origin[0] + size * x, origin[1] + size * y,
          origin[2] + size * z) for x, y, z in v]
    faces = [(0, 2, 1), (0, 3, 2), (4, 5, 6), (4, 6, 7),
             (0, 1, 5), (0, 5, 4), (3, 7, 6), (3, 6, 2),
             (0, 4, 7), (0, 7, 3), (1, 2, 6), (1, 6, 5)]
    if reverse:
        faces = [tuple(reversed(face)) for face in faces]
    return [[v[index] for index in face] for face in faces]


def _write_unit_cube_stl(path: Path) -> None:
    _write_ascii_stl(path, _unit_cube_triangles())


def test_missing_stl_is_input_error(tmp_path):
    r = _call(server.cubit_stl_to_vol, stl_path=str(tmp_path / "nope.stl"))
    assert r["status"] == "error" and r["kind"] == "input"


def test_bad_scheme_is_input_error(tmp_path):
    p = tmp_path / "x.stl"
    p.write_text("solid x\nendsolid x\n")
    r = _call(server.cubit_stl_to_vol, stl_path=str(p), scheme="prism")
    assert r["status"] == "error" and r["kind"] == "input"
    assert "hex" in r["error"] and "tet" in r["error"]


@pytest.mark.parametrize(
    ("kwargs", "message"),
    [
        ({"size": -1.0}, "size"),
        ({"size": float("nan")}, "size"),
        ({"closure_tolerance": -0.1}, "closure_tolerance"),
        ({"closure_tolerance": float("inf")}, "closure_tolerance"),
        ({"timeout_s": 0}, "timeout_s"),
        ({"timeout_s": float("nan")}, "timeout_s"),
    ],
)
def test_numeric_contract_fails_before_stl_inspection(tmp_path, kwargs,
                                                      message):
    p = tmp_path / "placeholder.stl"
    p.write_text("not inspected", encoding="ascii")
    r = _call(server.cubit_stl_to_vol, stl_path=str(p), **kwargs)
    assert r["status"] == "error" and r["kind"] == "input"
    assert message in r["error"]


def test_output_paths_must_not_alias_input_or_each_other(tmp_path):
    p = tmp_path / "cube.stl"
    _write_unit_cube_stl(p)

    same_as_input = _call(
        server.cubit_stl_to_vol, stl_path=str(p), out_vol=str(p))
    assert same_as_input["status"] == "error"
    assert "distinct paths" in same_as_input["error"]
    assert p.is_file(), "collision validation must not delete the source STL"

    shared = tmp_path / "shared.mesh"
    same_outputs = _call(
        server.cubit_stl_to_vol, stl_path=str(p),
        out_vol=str(shared), out_msh=str(shared))
    assert same_outputs["status"] == "error"
    assert "distinct paths" in same_outputs["error"]


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
    assert r["connected_components"] == 1
    assert abs(r["volume"] - 1.0) < 1e-12


def test_stl_volume_sums_disconnected_components_with_opposite_winding(
        tmp_path):
    p = tmp_path / "two_cubes.stl"
    triangles = [
        *_unit_cube_triangles(),
        *_unit_cube_triangles(origin=(2.0, 0.0, 0.0), reverse=True),
    ]
    _write_ascii_stl(p, triangles)

    r = inspect_stl(p)

    assert r["watertight"] is True
    assert r["winding_consistent"] is True
    assert r["connected_components"] == 2
    assert r["component_nesting_depths"] == [0, 0]
    assert r["volume"] == pytest.approx(2.0, rel=0.0, abs=1e-12)


def test_stl_volume_subtracts_a_nested_cavity_independent_of_winding(
        tmp_path):
    p = tmp_path / "hollow_cube.stl"
    triangles = [
        *_unit_cube_triangles(size=3.0),
        *_unit_cube_triangles(origin=(1.0, 1.0, 1.0), reverse=True),
    ]
    _write_ascii_stl(p, triangles)

    r = inspect_stl(p)

    assert r["watertight"] is True
    assert r["connected_components"] == 2
    assert r["component_nesting_depths"] == [0, 1]
    assert r["volume"] == pytest.approx(26.0, rel=0.0, abs=1e-12)


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


def test_mesh_total_volume_rejects_invalid_runtime_arguments(tmp_path):
    p = tmp_path / "placeholder.msh"
    p.write_text("placeholder", encoding="ascii")
    with pytest.raises(ValueError, match="quadrature"):
        mesh_total_volume(p, quadrature="")
    with pytest.raises(ValueError, match="timeout_s"):
        mesh_total_volume(p, timeout_s=float("inf"))


def test_headless_journal_rejects_invalid_timeout_before_environment_probe():
    with pytest.raises(ValueError, match="timeout_s"):
        server._cs.run_headless_journal(["reset"], timeout_s=float("nan"))


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


def test_relative_out_paths_are_absolutized(tmp_path, monkeypatch):
    """Relative out_vol/out_msh once resolved against themselves: the
    journal embeds the path while Cubit runs with cwd=vol.parent, so the
    .vol landed in a nested dir and the .msh writer failed with 'Cannot
    open file' (measured 2026-08-09 from a notebook kernel)."""
    p = tmp_path / "cube.stl"
    _write_unit_cube_stl(p)
    seen = {}

    def fake_headless(commands, **kwargs):
        seen["commands"] = list(commands)
        seen["kwargs"] = kwargs
        return {"status": "completed", "exit_code": 0}

    monkeypatch.setattr(server._cs, "run_headless_journal", fake_headless)
    monkeypatch.chdir(tmp_path)
    _call(server.cubit_stl_to_vol, stl_path=str(p),
          out_vol="work/cube.vol", out_msh="work/cube.msh")

    export_lines = [line for line in seen["commands"]
                    if line.startswith("export ")]
    assert len(export_lines) == 2
    root = server.PROJECT_ROOT.as_posix()
    for line in export_lines:
        quoted = line.split('"')[1]
        assert quoted.startswith(root), (quoted, root)
    assert Path(seen["kwargs"]["working_directory"]).is_absolute()


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
