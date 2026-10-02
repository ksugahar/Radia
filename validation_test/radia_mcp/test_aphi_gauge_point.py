"""Periodic A-Phi refuses a gauge that fixes no H1 dof, and accepts real ones."""

import pytest

ngsolve = pytest.importorskip("ngsolve")
occ = pytest.importorskip("netgen.occ")

from radia_mcp.radia_ngsolve.solve import _require_gauge_point


def _box(*, point="GND", isolated=False, face=None):
    box = occ.Box(occ.Pnt(0, 0, 0), occ.Pnt(1, 1, 1))
    if face:
        box.faces.Min(occ.X).name = face
    parts = [box]
    if point and isolated:
        gnd = occ.Vertex(occ.Pnt(3, 3, 3))
        gnd.name = point
        parts.append(gnd)
    elif point:
        box.vertices[0].name = point
    shape = occ.Glue(parts) if len(parts) > 1 else box
    with ngsolve.TaskManager():
        return ngsolve.Mesh(occ.OCCGeometry(shape).GenerateMesh(maxh=0.5))


def test_connected_gauge_point_is_accepted():
    _require_gauge_point(_box(), "", "GND")


def test_regular_expression_point_selector_is_accepted():
    _require_gauge_point(_box(point="REF"), "", "GND|REF")


@pytest.mark.parametrize("selector", ["REF", ""])
def test_point_selector_that_fixes_nothing_is_rejected(selector):
    with pytest.raises(ValueError, match="gauge is not fixed"):
        _require_gauge_point(_box(), "", selector)


def test_unconnected_gauge_vertex_is_rejected():
    mesh = _box(isolated=True)
    assert "GND" in set(mesh.GetBBBoundaries())
    with pytest.raises(ValueError, match="gauge is not fixed"):
        _require_gauge_point(mesh, "", "GND")


@pytest.mark.parametrize("selector", ["GND", ""])
def test_dirichlet_face_fixes_the_gauge_without_a_point(selector):
    _require_gauge_point(_box(point=None, face="outer"), "outer", selector)


def test_dirichlet_face_with_unconnected_point_is_accepted():
    _require_gauge_point(_box(isolated=True, face="outer"), "outer", "GND")
