"""Contract: periodic A-Phi refuses a gauge point that fixes no H1 vertex."""

import pytest

ngsolve = pytest.importorskip("ngsolve")
occ = pytest.importorskip("netgen.occ")

from radia_mcp.radia_ngsolve.solve import _require_gauge_point


def _mesh(*, isolated: bool):
    box = occ.Box(occ.Pnt(0, 0, 0), occ.Pnt(1, 1, 1))
    if isolated:
        gnd = occ.Vertex(occ.Pnt(3, 3, 3))
        gnd.name = "GND"
        shape = occ.Glue([box, gnd])
    else:
        box.vertices[0].name = "GND"
        shape = box
    return ngsolve.Mesh(occ.OCCGeometry(shape).GenerateMesh(maxh=0.5))


def test_connected_gauge_point_is_accepted():
    _require_gauge_point(_mesh(isolated=False), "", "GND")


def test_missing_gauge_tag_is_rejected():
    with pytest.raises(ValueError, match="not a point tag"):
        _require_gauge_point(_mesh(isolated=False), "", "REF")


def test_unconnected_gauge_vertex_is_rejected():
    mesh = _mesh(isolated=True)
    if "GND" not in set(mesh.GetBBBoundaries()):
        pytest.skip("this Netgen drops the free OCC vertex instead of tagging it")
    with pytest.raises(ValueError, match="fixes no H1 vertex"):
        _require_gauge_point(mesh, "", "GND")


def test_dirichlet_face_already_fixes_the_gauge():
    box = occ.Box(occ.Pnt(0, 0, 0), occ.Pnt(1, 1, 1))
    box.faces.Min(occ.X).name = "outer"
    gnd = occ.Vertex(occ.Pnt(3, 3, 3))
    gnd.name = "GND"
    mesh = ngsolve.Mesh(occ.OCCGeometry(occ.Glue([box, gnd])).GenerateMesh(maxh=0.5))
    if "GND" not in set(mesh.GetBBBoundaries()):
        pytest.skip("this Netgen drops the free OCC vertex instead of tagging it")
    _require_gauge_point(mesh, "outer", "GND")
