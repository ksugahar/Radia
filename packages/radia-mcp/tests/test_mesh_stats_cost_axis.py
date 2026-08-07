"""mesh_quality's cost-axis / anisotropy report.

The mesh-quality study (validation_test/radia_mcp/mesh_quality_study)
established that minSICN alone cannot rank meshes: the discriminating
axes are dof (not element count), the interior-node fraction, and
stretching. These lock the contract of the numbers that carry those.

The dof identities asserted here are structural, not conventions:
lowest-order H1 / HCurl / HDiv have exactly one dof per node / edge /
face, so the report must reproduce the counts an FE library would.
"""
from __future__ import annotations

import math
from pathlib import Path

import pytest

pytest.importorskip("gmsh")

from radia_mcp.gmsh.msh_inspect import mesh_quality


def _write_unit_cube_tet_msh(path: Path) -> None:
    """One cube split into 6 tets, written as .msh v4.1 by gmsh itself."""
    import gmsh

    gmsh.initialize(["-noconfig"])
    try:
        gmsh.option.setNumber("General.Terminal", 0)
        gmsh.model.add("cube")
        gmsh.model.occ.addBox(0, 0, 0, 1, 1, 1)
        gmsh.model.occ.synchronize()
        gmsh.option.setNumber("Mesh.MeshSizeMin", 1.0)
        gmsh.option.setNumber("Mesh.MeshSizeMax", 1.0)
        gmsh.model.mesh.generate(3)
        gmsh.write(str(path))
    finally:
        gmsh.finalize()


@pytest.fixture(scope="module")
def cube_msh(tmp_path_factory) -> Path:
    p = tmp_path_factory.mktemp("meshstats") / "cube.msh"
    _write_unit_cube_tet_msh(p)
    return p


def test_mesh_stats_present_and_self_consistent(cube_msh):
    q = mesh_quality(cube_msh)
    assert q["ran"], q
    stats = q["mesh_stats"]

    # dof_estimate is the whole point: nodes / edges / faces, i.e. the
    # dof count of the lowest-order H1 / HCurl / HDiv space.
    dof = stats["dof_estimate"]
    assert dof["h1_p1"] == stats["n_nodes"]
    assert dof["hcurl_lowest"] == stats["n_edges"]
    assert dof["hdiv_lowest"] == stats["n_faces_tri"] + stats["n_faces_quad"]
    assert dof["l2_p0"] == stats["n_elements_3d"]

    # a closed 3D mesh has boundary + interior nodes adding up
    assert (stats["n_interior_nodes"] + stats["n_boundary_nodes"]
            == stats["n_nodes"])
    assert 0.0 <= stats["interior_node_fraction"] <= 1.0
    assert math.isclose(
        stats["interior_node_fraction"],
        stats["n_interior_nodes"] / stats["n_nodes"], rel_tol=1e-12)

    # Euler check for a tet mesh of a ball-like (simply connected) solid:
    # V - E + F - C = 1
    euler = (stats["n_nodes"] - stats["n_edges"] + stats["n_faces_tri"]
             + stats["n_faces_quad"] - stats["n_elements_3d"])
    assert euler == 1, stats


def test_boundary_face_count_matches_face_bookkeeping(cube_msh):
    """Every interior face is shared by 2 cells, every boundary face by 1."""
    q = mesh_quality(cube_msh)
    stats = q["mesh_stats"]
    # tet mesh: 4 faces per cell
    elem_faces = 4 * stats["n_elements_3d"]
    unique = stats["n_faces_tri"]
    bnd = stats["n_boundary_faces"]
    assert elem_faces == 2 * unique - bnd, stats
    assert bnd > 0


def test_aspect_ratio_reported_and_at_least_one(cube_msh):
    q = mesh_quality(cube_msh)
    for bt in q["by_type"]:
        ar = bt["aspect_ratio"]
        # maxEdge/minEdge can never be below 1 for a non-degenerate element
        assert ar["min"] >= 1.0 - 1e-9, ar
        assert ar["min"] <= ar["mean"] <= ar["max"] + 1e-9
        assert ar["p95"] <= ar["max"] + 1e-9
        assert ar["n_nonfinite"] == 0
        iso = bt["min_isotropy"]
        assert 0.0 < iso["min"] <= 1.0 + 1e-9
        assert iso["min"] <= iso["mean"] + 1e-9


def test_stats_can_be_switched_off(cube_msh):
    q = mesh_quality(cube_msh, include_mesh_stats=False)
    assert q["ran"]
    assert "mesh_stats" not in q
    # the pre-existing quality contract is untouched either way
    assert q["metric"] == "minSICN"
    assert q["by_type"] and "min_quality" in q["by_type"][0]
    assert "aspect_ratio" not in q["by_type"][0]


