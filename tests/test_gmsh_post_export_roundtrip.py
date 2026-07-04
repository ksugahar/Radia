"""
GMSH .msh output round-trip validation (gmsh as an OUTPUT tool).

GmshPostExport (src/radia/gmsh_post_export.py) writes NGSolve mesh + field data
to GMSH .msh v4.1 with per-material physical groups and high-order elements.
These tests forge it as an *output* tool by re-opening every file we write with
the GMSH API itself and asserting the result is valid: node/element counts,
physical groups (by name), element TYPE codes (incl. high-order), and that the
field data lands as GMSH post-processing views. Pure NGSolve + gmsh (no Cubit).
"""
import os
import math
import pytest

ng = pytest.importorskip("ngsolve")
gmsh = pytest.importorskip("gmsh")
from ngsolve import CF, sqrt, x, y, z
from ngsolve.meshes import MakeStructured3DMesh

import sys
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
from radia.gmsh_post_export import GmshPostExport

_TMP = os.environ.get("TEMP", "C:\\temp")


def _reopen(path):
    """Open a .msh with the gmsh API; return (phys_named, vol_types, n_nodes, n_views)."""
    gmsh.initialize()
    try:
        gmsh.open(path)
        pg = gmsh.model.getPhysicalGroups()
        named = {gmsh.model.getPhysicalName(d, t) for d, t in pg}
        vtypes = list(gmsh.model.mesh.getElements(3)[0])
        n_nodes = len(gmsh.model.mesh.getNodes()[0])
        n_views = len(gmsh.view.getTags())
        return named, vtypes, n_nodes, n_views
    finally:
        gmsh.finalize()


def _assert_launch_companions(msh_path):
    geo_path = os.path.splitext(msh_path)[0] + ".geo"
    geo_opt_path = geo_path + ".opt"
    msh_opt_path = msh_path + ".opt"
    for path in (geo_path, geo_opt_path, msh_opt_path):
        assert os.path.exists(path), f"missing Gmsh launch companion: {path}"
    with open(geo_path, "r", encoding="utf-8") as fh:
        geo = fh.read()
    assert f'Merge "{os.path.basename(msh_path)}";' in geo
    assert "Mesh.NumSubEdges = 4;" in geo
    assert "General.RotationX = -68;" in geo
    return geo_path, geo_opt_path, msh_opt_path


def test_highorder_hex_roundtrip():
    """A curved (order-2) structured HEX mesh exports as Hex27 (gmsh type 12)
    and gmsh re-reads it with the right node count + 2 field views."""
    mesh = MakeStructured3DMesh(hexes=True, nx=3, ny=3, nz=3)
    mesh.Curve(2)
    post = GmshPostExport(mesh)
    post.add_vector_field("B", CF((y, -x, z)))
    post.add_scalar_field("absr", sqrt(x * x + y * y + z * z))
    out = os.path.join(_TMP, "rt_hex_test.msh")
    post.write(out)
    assert os.path.exists(out)
    _assert_launch_companions(out)
    named, vtypes, n_nodes, n_views = _reopen(out)
    assert 12 in vtypes, "expected Hex27 (type 12) in %s" % vtypes  # order-2 hex
    assert n_nodes == 343, "27 Hex27 share 343 nodes, got %d" % n_nodes
    assert n_views == 2, "expected 2 post views (B, absr), got %d" % n_views


def test_per_material_physical_groups():
    """A 2-material (sphere/air) curved tet mesh exports both materials as named
    GMSH physical groups, with high-order tets (Tet10, type 11)."""
    from netgen.csg import CSGeometry, OrthoBrick, Sphere, Pnt
    from ngsolve import Mesh
    geo = CSGeometry()
    geo.Add(Sphere(Pnt(0, 0, 0), 0.05).mat("sphere"))
    geo.Add(OrthoBrick(Pnt(-0.15, -0.15, -0.15), Pnt(0.15, 0.15, 0.15)).mat("air")
            - Sphere(Pnt(0, 0, 0), 0.05))
    mesh = Mesh(geo.GenerateMesh(maxh=0.06))
    mesh.Curve(2)
    post = GmshPostExport(mesh)
    post.add_vector_field("B", CF((y, -x, z)), material="sphere")
    post.add_scalar_field("phi", x + y + z)
    out = os.path.join(_TMP, "rt_2mat_test.msh")
    post.write(out)
    _assert_launch_companions(out)
    named, vtypes, n_nodes, n_views = _reopen(out)
    assert {"sphere", "air"} <= named, "physical groups %s missing sphere/air" % named
    assert 11 in vtypes or 4 in vtypes, "expected tet (type 11/4) in %s" % vtypes
    assert n_views == 2


def test_mesh_only_roundtrip():
    """write_mesh (no fields) still produces a gmsh-valid file with 0 views."""
    mesh = MakeStructured3DMesh(hexes=True, nx=2, ny=2, nz=2)
    post = GmshPostExport(mesh)
    out = os.path.join(_TMP, "rt_meshonly_test.msh")
    post.write_mesh(out)
    _assert_launch_companions(out)
    named, vtypes, n_nodes, n_views = _reopen(out)
    assert 5 in vtypes, "expected Hex8 (type 5) in %s" % vtypes  # order-1 hex
    assert n_nodes == 27
    assert n_views == 0


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
