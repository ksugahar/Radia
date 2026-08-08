"""Live regression for free Sculpt sidesets in the Netgen exporter.

This intentionally does *not* create mesh-based geometry.  The volume block
owns free hexes and the boundary owns direct/free quad faces, which exercises
the exporter path that geometry-only tests cannot reach.  Cubit is always run
in a disposable headless batch process.
"""
from __future__ import annotations

import os
from pathlib import Path

import pytest

pytest.importorskip("gmsh")
ng = pytest.importorskip("ngsolve")

from radia_mcp.cubit import session  # noqa: E402
from radia_mcp.cubit.server import _vol_surface_element_count  # noqa: E402
from radia_mcp.cubit.vol_inventory import summarize_netgen_vol_inventory  # noqa: E402
from radia_mcp.gmsh.msh_inspect import mesh_quality, mesh_total_volume  # noqa: E402
from cubit_mesh_export.check import check_consistency  # noqa: E402


pytestmark = pytest.mark.skipif(
    session.get_cubit_bin_dir() is None,
    reason="Coreform Cubit is not installed",
)


def _cube_stl_rows(x0: float, x1: float, name: str) -> list[str]:
    vertices = [
        (x0, 0, 0), (x1, 0, 0), (x1, 1, 0), (x0, 1, 0),
        (x0, 0, 1), (x1, 0, 1), (x1, 1, 1), (x0, 1, 1),
    ]
    faces = [
        (0, 2, 1), (0, 3, 2), (4, 5, 6), (4, 6, 7),
        (0, 1, 5), (0, 5, 4), (3, 7, 6), (3, 6, 2),
        (0, 4, 7), (0, 7, 3), (1, 2, 6), (1, 6, 5),
    ]
    rows = [f"solid {name}"]
    for face in faces:
        rows.extend(["  facet normal 0 0 0", "    outer loop"])
        rows.extend(
            f"      vertex {vertices[index][0]} {vertices[index][1]} "
            f"{vertices[index][2]}"
            for index in face
        )
        rows.extend(["    endloop", "  endfacet"])
    rows.append(f"endsolid {name}")
    return rows


def _write_cube_stl(path: Path) -> None:
    rows = _cube_stl_rows(0.0, 1.0, "cube")
    path.write_text("\n".join(rows) + "\n", encoding="ascii")


def _write_two_cube_stl(path: Path) -> None:
    rows = _cube_stl_rows(0.0, 1.0, "left")
    rows.extend(_cube_stl_rows(2.0, 3.0, "right"))
    path.write_text("\n".join(rows) + "\n", encoding="ascii")


def _write_two_material_box_stl(path: Path) -> None:
    path.write_text(
        "\n".join(_cube_stl_rows(0.0, 2.0, "two_material_box")) + "\n",
        encoding="ascii",
    )


def test_free_sculpt_sideset_reaches_vol_surfaceelements(tmp_path):
    stl = tmp_path / "cube.stl"
    vol = tmp_path / "cube.vol"
    msh = tmp_path / "cube.msh"
    _write_cube_stl(stl)

    commands = [
        f'import stl "{stl.as_posix()}" feature_angle 135 merge',
        "sculpt volume all processors 1 size 0.2 gen_sidesets 2",
        "block 1 add hex all",
        'block 1 name "solid"',
        f'export netgen "{vol.as_posix()}" order 2 overwrite',
        f'export gmsh "{msh.as_posix()}" dimension 3 order 1 overwrite',
    ]
    result = session.run_headless_journal(
        commands,
        timeout_s=900,
        working_directory=tmp_path,
        command_plugin_directory=os.environ.get("CUBIT_COMMAND_PLUGIN_DIR") or None,
    )

    assert result["status"] == "completed", result
    assert result["persistent_gui_started"] is False
    assert vol.is_file() and msh.is_file()
    vol_faces = _vol_surface_element_count(vol)
    quality = mesh_quality(msh)
    msh_faces = quality["mesh_stats"]["n_boundary_faces"]
    assert quality["total_negative"] == 0, quality
    assert vol_faces == msh_faces > 0
    curved_mesh = ng.Mesh(str(vol))
    assert curved_mesh.GetCurveOrder() == 2
    assert curved_mesh.ne == quality["mesh_stats"]["n_elements_3d"]
    volume = float(ng.Integrate(1, curved_mesh))
    outward_flux = float(ng.Integrate(
        ng.InnerProduct(
            ng.CF((ng.x, ng.y, ng.z)), ng.specialcf.normal(3)),
        curved_mesh, ng.BND))
    assert abs(outward_flux - 3.0 * volume) / (3.0 * volume) < 1e-12
    assert "added " in result["stdout_tail"]
    assert " free sideset faces under " in result["stdout_tail"]


def test_free_sculpt_sideset_recovers_each_material_domain(tmp_path):
    stl = tmp_path / "two_cubes.stl"
    vol = tmp_path / "two_cubes.vol"
    msh = tmp_path / "two_cubes.msh"
    _write_two_cube_stl(stl)

    commands = [
        f'import stl "{stl.as_posix()}" feature_angle 135 merge',
        "sculpt volume all processors 1 size 0.25 gen_sidesets 2",
        "delete block all",
        "block 1 add hex all with x_coord < 1.5",
        'block 1 name "left"',
        "block 2 add hex all with x_coord > 1.5",
        'block 2 name "right"',
        f'export netgen "{vol.as_posix()}" order 2 overwrite',
        f'export gmsh "{msh.as_posix()}" dimension 3 order 1 overwrite',
    ]
    result = session.run_headless_journal(
        commands,
        timeout_s=900,
        working_directory=tmp_path,
        command_plugin_directory=os.environ.get("CUBIT_COMMAND_PLUGIN_DIR") or None,
    )

    assert result["status"] == "completed", result
    assert result["persistent_gui_started"] is False
    inventory = summarize_netgen_vol_inventory(
        vol.read_text(encoding="utf-8"), source=str(vol)
    )
    ownership = inventory["boundary_domain_ownership"]
    assert ownership["passed"], ownership
    assert ownership["unreferenced_volume_domains"] == []
    assert set(ownership["domain_surface_incidence_counts"]) == {1, 2}
    assert ownership["domain_surface_incidence_counts"][1] > 0
    assert (
        ownership["domain_surface_incidence_counts"][1]
        == ownership["domain_surface_incidence_counts"][2]
    )
    assert ownership["surface_domain_pair_counts"] == {
        f"1->0": ownership["domain_surface_incidence_counts"][1],
        f"2->0": ownership["domain_surface_incidence_counts"][2],
    }

    mesh = ng.Mesh(str(vol))
    assert set(mesh.GetMaterials()) == {"left", "right"}
    left_volume = float(ng.Integrate(1, mesh, definedon=mesh.Materials("left")))
    right_volume = float(ng.Integrate(1, mesh, definedon=mesh.Materials("right")))
    assert left_volume > 0.0
    assert left_volume == pytest.approx(right_volume, rel=1e-6, abs=1e-9)
    gmsh_volume = float(mesh_total_volume(msh)["total_volume"])
    assert left_volume + right_volume == pytest.approx(
        gmsh_volume, rel=5e-6, abs=1e-9
    )
    total_flux = float(ng.Integrate(
        ng.InnerProduct(ng.CF((ng.x, ng.y, ng.z)), ng.specialcf.normal(3)),
        mesh,
        ng.BND,
    ))
    assert abs(total_flux - 3.0 * (left_volume + right_volume)) < 1e-11

    report = check_consistency(vol)
    assert report["passed"], report["warnings"]
    assert report["boundary_domain_ownership"]["passed"]
    ng_area = sum(row["ng_area"] for row in report["boundaries"])
    cad_area = sum(row.get("cad_area", 0.0) for row in report["boundaries"])
    assert ng_area > 0.0
    area_error_pct = abs(ng_area - cad_area) / cad_area * 100.0
    assert area_error_pct < 0.05
    assert cad_area == pytest.approx(ng_area, rel=5e-4, abs=1e-12)
    assert "domain-aware synthetic descriptors" in result["stdout_tail"]


def test_free_sculpt_skin_preserves_internal_material_interface_once(tmp_path):
    stl = tmp_path / "two_material_box.stl"
    vol = tmp_path / "two_material_box.vol"
    msh = tmp_path / "two_material_box.msh"
    _write_two_material_box_stl(stl)

    commands = [
        f'import stl "{stl.as_posix()}" feature_angle 135 merge',
        "sculpt volume all processors 1 size 0.25 gen_sidesets 2",
        "delete block all",
        "block 1 add hex all with x_coord < 1.0",
        'block 1 name "left"',
        "block 2 add hex all with x_coord > 1.0",
        'block 2 name "right"',
        "skin block 1 make sideset 3",
        'sideset 3 name "left_right_interface"',
        f'export netgen "{vol.as_posix()}" order 2 overwrite',
        f'export gmsh "{msh.as_posix()}" dimension 3 order 1 overwrite',
    ]
    result = session.run_headless_journal(
        commands,
        timeout_s=900,
        working_directory=tmp_path,
        command_plugin_directory=os.environ.get("CUBIT_COMMAND_PLUGIN_DIR") or None,
    )

    assert result["status"] == "completed", result
    assert result["persistent_gui_started"] is False
    inventory = summarize_netgen_vol_inventory(
        vol.read_text(encoding="utf-8"), source=str(vol)
    )
    ownership = inventory["boundary_domain_ownership"]
    assert ownership["passed"], ownership
    assert ownership["internal_interface_element_count"] > 0
    assert ownership["internal_interface_domain_pair_counts"] == {
        "1<->2": ownership["internal_interface_element_count"]
    }
    assert (
        ownership["exterior_surface_element_count"]
        + ownership["internal_interface_element_count"]
        == inventory["surface_elements"]
    )
    assert ownership["duplicate_surface_connectivity_rows"] == []

    mesh = ng.Mesh(str(vol))
    interface = mesh.Boundaries("left_right_interface")
    exterior = mesh.Boundaries("sideset_1")
    interface_area = float(
        ng.Integrate(1, mesh, ng.BND, definedon=interface, order=16)
    )
    assert interface_area > 0.0
    volume = float(ng.Integrate(1, mesh))
    exterior_flux = float(ng.Integrate(
        ng.InnerProduct(ng.CF((ng.x, ng.y, ng.z)), ng.specialcf.normal(3)),
        mesh,
        ng.BND,
        definedon=exterior,
    ))
    assert exterior_flux == pytest.approx(3.0 * volume, rel=1e-12, abs=1e-11)

    report = check_consistency(vol, threshold=0.05)
    assert report["passed"], report["warnings"]
    by_name = {row["name"]: row for row in report["boundaries"]}
    assert by_name["left_right_interface"]["cad_area"] == pytest.approx(
        interface_area, rel=1e-12, abs=1e-12
    )
    assert abs(by_name["left_right_interface"]["error_pct"]) < 1e-9
    assert "added " in result["stdout_tail"]
    assert "domain-aware synthetic descriptors" in result["stdout_tail"]
