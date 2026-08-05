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
import numpy as np
import pytest

ng = pytest.importorskip("ngsolve")
gmsh = pytest.importorskip("gmsh")
from ngsolve import CF, sqrt, x, y, z
from ngsolve.meshes import MakeStructured3DMesh

import sys
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
from radia.gmsh_post_export import (
    GmshPostExport,
    _get_gmsh_ref_points,
    export_element_activation_animation,
    export_nodal_deformation_animation,
)


def test_gmsh_reference_points_use_element_dimension_stride():
    trig_type, trig = _get_gmsh_ref_points("TRIG", 2)
    quad_type, quad = _get_gmsh_ref_points("QUAD", 2)

    assert trig_type is not None and len(trig) == 6
    assert quad_type is not None and len(quad) == 9
    assert all(len(point) == 3 for point in trig + quad)

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


def _assert_positive_jacobians(path, expect_vol=None, rel=1e-9):
    """Gate: no negative getJacobians point on any 3D element; optionally
    check the integrated volume against the exact value (repo policy:
    negative determinants = inverted node ordering)."""
    gmsh.initialize()
    try:
        gmsh.open(path)
        total_neg = 0
        vol = 0.0
        etypes, _, _ = gmsh.model.mesh.getElements(3)
        assert etypes, "no 3D elements in %s" % path
        for et in etypes:
            local, weights = gmsh.model.mesh.getIntegrationPoints(
                int(et), "Gauss2")
            _, det, _ = gmsh.model.mesh.getJacobians(int(et), local)
            det = np.asarray(det, dtype=float)
            w = np.asarray(weights, dtype=float)
            total_neg += int((det <= 0.0).sum())
            vol += float((det.reshape(-1, len(w)) * w).sum())
    finally:
        gmsh.finalize()
    assert total_neg == 0, f"{total_neg} negative Jacobian points in {path}"
    if expect_vol is not None:
        assert vol == pytest.approx(expect_vol, rel=rel)


def _shear(x, y, z):
    # det J = 1 everywhere, so the mapped cube keeps volume exactly 1,
    # while every element becomes asymmetric (catches orientation and
    # high-order node placement bugs that symmetric cubes can hide).
    return (x + 0.3 * y * y, y, z + 0.2 * x)


def test_volume_element_orientation_tet_hex_prism():
    """TET/HEX/PRISM exports must have strictly positive Jacobians and
    exact volumes (locks the 2026-08 orientation fixes: TET corner
    permutation [3,0,1,2], HEX y/z-swapped and PRISM barycentric-cycled
    high-order reference transforms)."""
    from netgen.occ import OCCGeometry, Sphere
    from netgen.occ import Pnt as OccPnt
    from ngsolve import Mesh

    geo = OCCGeometry(Sphere(OccPnt(0, 0, 0), 1.0))
    mesh = Mesh(geo.GenerateMesh(maxh=0.5))
    mesh.Curve(2)
    out = os.path.join(_TMP, "rt_orient_tet.msh")
    GmshPostExport(mesh).write_mesh(out)
    _assert_positive_jacobians(out, expect_vol=4.0 / 3.0 * math.pi, rel=5e-3)

    for tag, kwargs in (("hex", dict(hexes=True)),
                        ("prism", dict(hexes=False, prism=True))):
        mesh = MakeStructured3DMesh(nx=2, ny=2, nz=2, mapping=_shear,
                                    **kwargs)
        mesh.Curve(2)
        out = os.path.join(_TMP, f"rt_orient_{tag}.msh")
        GmshPostExport(mesh).write_mesh(out)
        _assert_positive_jacobians(out, expect_vol=1.0)


def test_volume_element_orientation_pyramid():
    """An asymmetric hand-built PYR5/PYR14 export has positive Jacobians
    and the exact cone volume (locks the identity corner permutation and
    the frustum-bijection reference transform)."""
    from netgen.csg import Pnt
    from netgen.meshing import Element2D, Element3D, FaceDescriptor
    from netgen.meshing import Mesh as NgMesh
    from netgen.meshing import MeshPoint
    from ngsolve import Mesh

    m = NgMesh(dim=3)
    pts = [(0, 0, 0), (1.2, 0, 0), (1.1, 1.3, 0), (0, 1, 0),
           (0.3, 0.2, 1.0)]
    pids = [m.Add(MeshPoint(Pnt(*pt))) for pt in pts]
    m.SetMaterial(1, "pyr")
    m.Add(FaceDescriptor(surfnr=1, domin=1, bc=1))
    m.Add(Element3D(1, pids))
    m.Add(Element2D(1, [pids[3], pids[2], pids[1], pids[0]]))
    for a, b in ((0, 1), (1, 2), (2, 3), (3, 0)):
        m.Add(Element2D(1, [pids[a], pids[b], pids[4]]))
    base_area = 1.33  # shoelace of the planar z=0 quad
    for order, tag in ((1, "pyr1"), (2, "pyr2")):
        mesh = Mesh(m)
        if order >= 2:
            mesh.Curve(order)
        out = os.path.join(_TMP, f"rt_orient_{tag}.msh")
        GmshPostExport(mesh).write_mesh(out)
        _assert_positive_jacobians(out, expect_vol=base_area / 3.0, rel=1e-6)


def _probe_view(path, points):
    gmsh.initialize()
    try:
        gmsh.open(path)
        tag = gmsh.view.getTags()[0]
        out = []
        for pt in points:
            vals, dist = gmsh.view.probe(tag, *pt)
            assert dist == 0.0, f"probe point {pt} outside data"
            out.append(np.asarray(vals, dtype=float))
    finally:
        gmsh.finalize()
    return out


def test_curved_nodal_cf_field_probes_exactly():
    """CF fields are evaluated at EVERY emitted node of a curved mesh.

    Locks the 2026-08 vertex-only NodeData bug: with mid-edge values
    missing, gmsh's quadratic interpolation returned garbage (vertex
    shape functions are negative at edge midpoints). A linear vector
    field must now probe back exactly anywhere in the mesh."""
    from netgen.occ import OCCGeometry, Sphere
    from netgen.occ import Pnt as OccPnt
    from ngsolve import Mesh

    geo = OCCGeometry(Sphere(OccPnt(0, 0, 0), 1.0))
    mesh = Mesh(geo.GenerateMesh(maxh=0.4))
    mesh.Curve(2)
    post = GmshPostExport(mesh)
    post.add_vector_field("B", CF((y, -x, z)))
    out = os.path.join(_TMP, "rt_curved_field_exact.msh")
    post.write(out)

    points = [(0.5, 0.0, 0.0), (0.1, 0.2, 0.3), (-0.3, 0.4, -0.2)]
    for pt, vals in zip(points, _probe_view(out, points)):
        np.testing.assert_allclose(vals[:3], (pt[1], -pt[0], pt[2]),
                                   atol=1e-9)


@pytest.mark.parametrize("kwargs", (dict(hexes=True),
                                    dict(hexes=False),
                                    dict(hexes=False, prism=True)))
def test_vertex_array_field_expands_exactly(kwargs):
    """Vertex-length arrays are embedded exactly onto order-2 nodes.

    On straight-edged Curve(2) meshes the extra Lagrange nodes sit at
    corner means, so a linear field given per-vertex must probe back
    exactly -- this also locks the hex27 face/center and prism18 face
    orderings of the P1->P2 fill tables."""
    mesh = MakeStructured3DMesh(nx=2, ny=2, nz=2, **kwargs)
    mesh.Curve(2)
    data = np.array([1.0 + 2.0 * v.point[0] - v.point[1]
                     + 0.5 * v.point[2] for v in mesh.vertices])
    post = GmshPostExport(mesh)
    post.add_scalar_field("phi", data)
    tag = "".join(k for k in kwargs)
    out = os.path.join(_TMP, f"rt_vertex_expand_{tag}.msh")
    post.write(out)

    points = [(0.51, 0.24, 0.77), (0.13, 0.88, 0.31)]
    for pt, vals in zip(points, _probe_view(out, points)):
        expect = 1.0 + 2.0 * pt[0] - pt[1] + 0.5 * pt[2]
        np.testing.assert_allclose(vals[0], expect, atol=1e-9)


def test_export_passes_radia_mcp_verify_gate():
    """Shipping gate: a fresh curved export passes the full radia-mcp
    verification (structural + NaN/Inf + Jacobian gates on the .msh
    plus the deep launch check on the companion .geo)."""
    verify_mod = pytest.importorskip("radia_mcp.gmsh.verify")

    mesh = MakeStructured3DMesh(hexes=True, nx=2, ny=2, nz=2)
    mesh.Curve(2)
    post = GmshPostExport(mesh)
    post.add_scalar_field("phi", x + y + z)
    out = os.path.join(_TMP, "rt_mcp_verify_gate.msh")
    post.write(out)

    result = verify_mod.verify_artifact(out, check_jacobians=True)
    assert result["ok"], result
    assert f"msh:{os.path.basename(out)}" in result["passed"]
    geo_gates = [g for g in result["passed"] if g.startswith("geo:")]
    assert geo_gates, "companion .geo gate did not run"


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


def test_element_activation_animation_roundtrip():
    mesh = MakeStructured3DMesh(hexes=True, nx=2, ny=1, nz=1)
    out = os.path.join(_TMP, "rt_element_activation_animation.msh")
    artifact = export_element_activation_animation(
        mesh, [[True, False], [True, True], [True, True]], out)
    _assert_launch_companions(out)
    assert artifact["active_counts"] == [1, 2, 2]

    gmsh.initialize()
    try:
        gmsh.open(out)
        tags = gmsh.view.getTags()
        assert len(tags) == 1
        for step, expected in enumerate((1, 2, 2)):
            kind, element_tags, values, time, ncomp = (
                gmsh.view.getModelData(tags[0], step))
            assert kind == "ElementData"
            assert len(element_tags) == expected
            assert len(values) == expected
            assert time == pytest.approx(float(step))
            assert ncomp == 1
    finally:
        gmsh.finalize()


def test_element_activation_animation_rejects_removal():
    mesh = MakeStructured3DMesh(hexes=True, nx=2, ny=1, nz=1)
    with pytest.raises(ValueError, match="removes active elements"):
        export_element_activation_animation(
            mesh, [[True, True], [True, False]],
            os.path.join(_TMP, "rt_nonmonotone_activation.msh"))


@pytest.mark.parametrize("invalid", ([True, 0.5], [True, np.nan]))
def test_element_activation_animation_rejects_nonbinary_masks(invalid):
    mesh = MakeStructured3DMesh(hexes=True, nx=2, ny=1, nz=1)
    with pytest.raises(ValueError, match="boolean or binary"):
        export_element_activation_animation(
            mesh, [invalid],
            os.path.join(_TMP, "rt_invalid_activation.msh"))


def test_nodal_deformation_animation_roundtrip():
    mesh = MakeStructured3DMesh(hexes=True, nx=1, ny=1, nz=1)
    nodes = np.asarray([vertex.point for vertex in mesh.vertices], dtype=float)
    displacement = np.column_stack((
        0.01 * nodes[:, 0], -0.02 * nodes[:, 1], 0.03 * nodes[:, 2]))
    out = os.path.join(_TMP, "rt_nodal_deformation_animation.msh")
    artifact = export_nodal_deformation_animation(
        mesh, [np.zeros_like(displacement), displacement], out)
    geo_path, _, _ = _assert_launch_companions(out)
    assert artifact["steps"] == 2
    assert artifact["maximum_displacement_m"][1] > 0.0
    with open(geo_path, "r", encoding="utf-8") as stream:
        assert "View[0].VectorType = 5;" in stream.read()

    gmsh.initialize()
    try:
        gmsh.open(out)
        tags = gmsh.view.getTags()
        assert len(tags) == 1
        kind, node_tags, values, time_value, ncomp = (
            gmsh.view.getModelData(tags[0], 1))
        assert kind == "NodeData"
        assert ncomp == 3
        assert time_value == pytest.approx(1.0)
        assert len(node_tags) == len(nodes)
        np.testing.assert_allclose(np.asarray(values), displacement)
    finally:
        gmsh.finalize()


def test_curved_nodal_deformation_samples_ngsolve_callable():
    mesh = MakeStructured3DMesh(hexes=True, nx=1, ny=1, nz=1)
    mesh.Curve(2)
    displacement = CF((0.01 * x, -0.02 * y, 0.03 * z))
    out = os.path.join(_TMP, "rt_curved_nodal_deformation.msh")
    export_nodal_deformation_animation(mesh, [displacement], out)

    gmsh.initialize()
    try:
        gmsh.open(out)
        tags = gmsh.view.getTags()
        assert len(tags) == 1
        kind, node_tags, values, _, ncomp = gmsh.view.getModelData(tags[0], 0)
        assert kind == "NodeData" and ncomp == 3
        assert len(node_tags) > mesh.nv
        assert np.all(np.isfinite(np.asarray(values)))
    finally:
        gmsh.finalize()


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
