"""The rotational HDiv seam tie must not depend on vertex numbering.

Guards: on one straight HEX sector cell (fold 12) with three vertex numberings of the same geometry, the reduced
``Compress(Periodic(HDiv))`` space reproduces the seam-varying periodic field f = (-y, x, 0) (1e5/0.04) (1 + 50 z) to
round-off at BDM1 and BDM2. Its exact squared L2 norm (32914.709526) is integrated on the mesh, independently of the
reduced space. A numbering whose seam faces order their vertices unlike their image is refused before NGSolve ties it.
``renumber_periodic_seam`` repairs such a mesh and moves slave vertices only. With the stock tie, the swapped numbering
loses 1.5e-2 (p1) and 8.5e-3 (p2) of the norm, because NGSolve ties its slave face to a non-seam face.
Assumes: the field is in both spaces on this straight cell; exact representability is the oracle, so 1e-12 is a
round-off bound, not an accuracy criterion. The check is conservative: it also refuses the shuffled numbering, whose
stock tie happens to be right. Only PERIODIC identifications count, an empty seam refuses, and a CAD-bound mesh
is not renumbered.
"""

import math

import numpy as np
import pytest

ng = pytest.importorskip("ngsolve")

SEAM = ("periodic_min", "periodic_max")
NUMBERINGS = dict(sector=None, shuffled=[5, 0, 1, 7, 2, 6, 4, 3], swapped=[2, 6, 7, 5, 1, 3, 0, 4])
EXACT = 32914.709526
ROUND_OFF = 1e-12


def _cell(numbering, fold=12, inner=0.03, outer=0.05, half=0.004, identified=True):
    """One straight annular HEX cell; `numbering` holds the vertex number of each (station, radius, height)."""
    from netgen.meshing import Element2D, Element3D, FaceDescriptor, IdentificationType, Mesh, MeshPoint, Pnt
    mesh = Mesh(dim=3)
    mesh.SetMaterial(1, "yoke")
    skin = mesh.Add(FaceDescriptor(surfnr=1, domin=1, domout=0, bc=1))
    low = mesh.Add(FaceDescriptor(surfnr=2, domin=1, domout=0, bc=2))
    high = mesh.Add(FaceDescriptor(surfnr=3, domin=1, domout=0, bc=3))
    for index, name in enumerate(("skin",) + SEAM):
        mesh.SetBCName(index, name)
    keys = list(np.ndindex(2, 2, 2))
    point = {}
    step = 2.0*math.pi/fold
    for it, ir, iz in (keys if numbering is None else [keys[k] for k in np.argsort(numbering)]):
        radius = (inner, outer)[ir]
        point[it, ir, iz] = mesh.Add(MeshPoint(Pnt(radius*math.cos(it*step), radius*math.sin(it*step), (-half, half)[iz])))
    corners = [point[0, 0, 0], point[0, 1, 0], point[1, 1, 0], point[1, 0, 0], point[0, 0, 1], point[0, 1, 1], point[1, 1, 1], point[1, 0, 1]]
    mesh.Add(Element3D(1, corners))
    for quad in ((0, 3, 2, 1), (4, 5, 6, 7), (1, 2, 6, 5), (3, 0, 4, 7)):
        mesh.Add(Element2D(skin, [corners[k] for k in quad]))
    mesh.Add(Element2D(low, [corners[k] for k in (0, 1, 5, 4)]))
    mesh.Add(Element2D(high, [corners[k] for k in (2, 3, 7, 6)]))
    for ir, iz in np.ndindex(2, 2) if identified else ():
        mesh.AddPointIdentification(point[0, ir, iz], point[1, ir, iz], identnr=1, type=IdentificationType.PERIODIC)
    return ng.Mesh(mesh)


def _field():
    return ng.CF((-ng.y, ng.x, 0))*(1e5/0.04)*(1 + 50*ng.z)


def _deficit(mesh, order):
    from radia.vim._solve import _hdiv_space_with_image_constraints
    fes = _hdiv_space_with_image_constraints(mesh, order, [], SEAM)[0]
    u, v = fes.TnT()
    rule = {ng.HEX: ng.IntegrationRule(ng.HEX, 12)}
    mass = ng.BilinearForm(u*v*ng.dx(intrules=rule)).Assemble()
    load = ng.LinearForm(_field()*v*ng.dx(intrules=rule)).Assemble()
    r = np.array(load.vec.FV().NumPy())
    exact = float(ng.Integrate(ng.InnerProduct(_field(), _field()), mesh, order=12))
    return exact, (exact - float(r @ np.linalg.solve(np.array(mass.mat.ToDense()), r)))/exact


@pytest.mark.parametrize("name", ["shuffled", "swapped"])
def test_a_misoriented_seam_numbering_refuses(name):
    from radia.vim._solve import _hdiv_space_with_image_constraints
    with pytest.raises(ValueError, match="order their vertex numbers"):
        _hdiv_space_with_image_constraints(_cell(NUMBERINGS[name]), 2, [], SEAM)


@pytest.mark.parametrize("order", [1, 2])
@pytest.mark.parametrize("name", list(NUMBERINGS))
def test_the_renumbered_seam_reproduces_a_periodic_field(name, order):
    from radia.vim import renumber_periodic_seam
    raw = _cell(NUMBERINGS[name])
    fixed = renumber_periodic_seam(raw, SEAM)
    slaves = {int(getattr(b, "nr", b)) for _a, b, _n in raw.ngmesh.GetIdentifications()}
    moved = [k + 1 for k, (p, q) in enumerate(zip(raw.ngmesh.Coordinates(), fixed.ngmesh.Coordinates())) if not np.array_equal(p, q)]
    assert set(moved) <= slaves, f"{name}: non-slave vertices moved: {moved} slaves {sorted(slaves)}"
    exact, deficit = _deficit(fixed, order)
    assert abs(exact - EXACT) <= 1e-6*EXACT, f"exact norm {exact!r} expected {EXACT!r}"
    assert abs(deficit) <= ROUND_OFF, f"{name} order={order} relative projection deficit={deficit:.3e} limit={ROUND_OFF:.0e}"


def test_the_seam_contract_refuses_other_identification_types_empty_seams_and_cad_bound_meshes():
    from netgen.meshing import IdentificationType
    from netgen.occ import Box, OCCGeometry
    from radia.vim import check_periodic_seam_orientation, renumber_periodic_seam
    closing = _cell(None)
    pairs = closing.ngmesh.GetIdentifications()
    closing = _cell(None, identified=False)
    for a, b, _n in pairs:
        closing.ngmesh.AddPointIdentification(a, b, identnr=1, type=IdentificationType.CLOSESURFACES)
    closing = ng.Mesh(closing.ngmesh)
    with pytest.raises(ValueError, match="no PERIODIC partner"):
        check_periodic_seam_orientation(closing, SEAM)
    with pytest.raises(ValueError, match="have no faces"):
        check_periodic_seam_orientation(_cell(None), ("skin_a", "skin_b"))
    cad = ng.Mesh(OCCGeometry(Box((0, 0, 0), (1, 1, 1))).GenerateMesh(maxh=0.5))
    with pytest.raises(ValueError, match="OCCGeometry"):
        renumber_periodic_seam(cad, SEAM)
