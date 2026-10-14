"""Asking the mesh beats rebuilding what it already knows.

`boundaries_touching_materials` answers "which of these boundaries has a face
on an element of one of these materials?" with region algebra. The obvious
alternative -- walk every volume element, then every face of every element,
into a dictionary of sets -- is what the solver used to do, and it cost 2.05 s
on a 274k-element mesh, about a quarter of the solve it sat in.

So these tests pin the answer against that exhaustive definition rather than
against a stored expectation: if the region form ever stops meaning the same
thing, the brute force here disagrees with it.
"""
from __future__ import annotations

import pytest

pytestmark = pytest.mark.usefixtures("ngsolve_taskmanager")


def _brute_force(mesh, names, materials):
    """The exhaustive definition, written out. Correct, and slow by design."""
    import ngsolve as ng

    wanted = {str(n) for n in names}
    material_set = {str(m) for m in materials}
    face_materials = {}
    for element in mesh.Elements(ng.VOL):
        for face in element.faces:
            face_materials.setdefault(face.nr, set()).add(element.mat)
    touching = set()
    for element in mesh.Elements(ng.BND):
        if element.mat not in wanted:
            continue
        adjacent = set()
        for face in element.faces:
            adjacent.update(face_materials.get(face.nr, ()))
        if adjacent & material_set:
            touching.add(element.mat)
    return touching


@pytest.fixture(scope="module")
def two_material_mesh():
    """A box of one material stacked on a box of another, all faces named."""
    ng = pytest.importorskip("ngsolve")
    from netgen.occ import Z, Box, Glue, OCCGeometry, Pnt

    lower = Box(Pnt(0, 0, 0), Pnt(1, 1, 1))
    lower.mat("lower")
    lower.faces.Min(Z).name = "bottom"
    upper = Box(Pnt(0, 0, 1), Pnt(1, 1, 2))
    upper.mat("upper")
    upper.faces.Max(Z).name = "top"
    shape = Glue([lower, upper])
    for face in shape.faces:
        if not face.name:
            face.name = "side"
    return ng.Mesh(OCCGeometry(shape).GenerateMesh(maxh=0.35))


def test_it_agrees_with_walking_every_face(two_material_mesh):
    """The region answer is the exhaustive answer, on a real mesh."""
    from radia.kelvin_solver import boundaries_touching_materials

    mesh = two_material_mesh
    names = sorted(set(mesh.GetBoundaries()))
    assert names, "the fixture lost its boundary names"
    for materials in (("lower",), ("upper",), ("lower", "upper")):
        fast = boundaries_touching_materials(mesh, names, materials)
        slow = _brute_force(mesh, names, materials)
        assert fast == slow, (materials, sorted(fast), sorted(slow))


def test_the_geometry_it_reports_is_the_geometry_that_is_there(
        two_material_mesh):
    """Not merely self-consistent: the answers match the stacked boxes."""
    from radia.kelvin_solver import boundaries_touching_materials

    mesh = two_material_mesh
    assert "bottom" in boundaries_touching_materials(
        mesh, ("bottom", "top"), ("lower",))
    assert "top" not in boundaries_touching_materials(
        mesh, ("bottom", "top"), ("lower",))
    assert "top" in boundaries_touching_materials(
        mesh, ("bottom", "top"), ("upper",))
    assert "bottom" not in boundaries_touching_materials(
        mesh, ("bottom", "top"), ("upper",))
    # The sides run the full height, so they touch both.
    assert "side" in boundaries_touching_materials(mesh, ("side",), ("lower",))
    assert "side" in boundaries_touching_materials(mesh, ("side",), ("upper",))


def test_unknown_names_are_ignored_not_raised(two_material_mesh):
    """A scan over the mesh's own elements would simply never match them."""
    from radia.kelvin_solver import boundaries_touching_materials

    mesh = two_material_mesh
    got = boundaries_touching_materials(
        mesh, ("bottom", "a-boundary-that-is-not-here"), ("lower",))
    assert got == {"bottom"}


def test_an_unknown_material_selects_nothing(two_material_mesh):
    """An empty material selection cannot silently mean 'everything'."""
    from radia.kelvin_solver import boundaries_touching_materials

    mesh = two_material_mesh
    names = sorted(set(mesh.GetBoundaries()))
    assert boundaries_touching_materials(mesh, names, ("not-a-material",)) == set()
    assert boundaries_touching_materials(mesh, names, ()) == set()
