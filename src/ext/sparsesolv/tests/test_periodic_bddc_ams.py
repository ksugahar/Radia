"""Guards: periodic HCurl edge numbering in wirebasket AMS. Assumes: two translated tetrahedra."""

import numpy as np
import pytest

ng = pytest.importorskip("ngsolve")


def test_periodic_wirebasket_ams_matches_sparsecholesky():
    """The paired face has three slave edge DOFs and permuted vertex numbers."""
    from netgen.csg import Pnt
    from netgen.meshing import (
        Element2D, Element3D, FaceDescriptor, IdentificationType, Mesh, MeshPoint,
    )
    from ngsolve.krylovspace import CGSolver
    from radia import sparsesolv_ngsolve as ssn

    native = Mesh(dim=3)
    points = [
        native.Add(MeshPoint(Pnt(*xyz)))
        for xyz in ((0, 0, 0), (1, 0, 0), (0, 1, 0), (0, 0, 1),
                    (3, 1, 0), (3, 0, 0), (4, 0, 0), (3, 0, 1))
    ]
    for material, indices, name in (
        (1, (0, 1, 2, 3), "air"), (2, (5, 6, 4, 7), "kelvin"),
    ):
        native.SetMaterial(material, name)
        native.Add(Element3D(material, [points[i] for i in indices]))
        boundary = native.Add(
            FaceDescriptor(surfnr=material, domin=material, domout=0, bc=material)
        )
        for face in ((0, 2, 1), (0, 1, 3), (0, 3, 2), (1, 2, 3)):
            native.Add(Element2D(boundary, [points[indices[i]] for i in face]))
    for master, slave in ((0, 5), (1, 6), (2, 4)):
        native.AddPointIdentification(
            points[master], points[slave], 1, IdentificationType.PERIODIC
        )

    mesh = ng.Mesh(native)
    fes = ng.Periodic(ng.HCurl(mesh, order=2, nograds=True))
    free = fes.FreeDofs()
    assert sum(free[edge] for edge in range(mesh.nedge)) == mesh.nedge - 3
    for dof in range(mesh.nedge, fes.ndof):
        if fes.CouplingType(dof) == ng.COUPLING_TYPE.WIREBASKET_DOF:
            fes.SetCouplingType(dof, ng.COUPLING_TYPE.INTERFACE_DOF)

    u, v = fes.TnT()
    r2 = (ng.x - 3)**2 + ng.y**2 + ng.z**2
    curl_form = ng.curl(u) * ng.curl(v)
    mass_form = 1e-6 * u * v
    matrix = ng.BilinearForm(fes, symmetric=True)
    matrix += (curl_form + mass_form) * ng.dx("air")
    matrix += (r2 * curl_form + mass_form / (r2 + 1e-20)) * ng.dx(
        "kelvin", bonus_intorder=4
    )
    rhs = ng.LinearForm((ng.curl(v) * ng.CF((1, 2, 3))) * ng.dx("air"))
    pre = ng.Preconditioner(
        matrix, "bddc", coarsetype="sparsesolv_ams", coarseflags={"cycles": 3}
    )
    with ng.TaskManager():
        matrix.Assemble()
        rhs.Assemble()
    direct = ng.GridFunction(fes)
    direct.vec.data = matrix.mat.Inverse(free, inverse="sparsecholesky") * rhs.vec
    ams = ng.GridFunction(fes)
    cg = CGSolver(matrix.mat, pre.mat, tol=1e-12, maxiter=300, printrates=False)
    ams.vec.data = cg * rhs.vec
    residual = rhs.vec.CreateVector()
    residual.data = rhs.vec - matrix.mat * ams.vec
    mask = np.array(free, dtype=bool)
    assert np.linalg.norm(residual.FV().NumPy()[mask]) <= (
        1e-6 * np.linalg.norm(rhs.vec.FV().NumPy()[mask])
    )
    point = mesh(0.2, 0.2, 0.2)
    np.testing.assert_allclose(ng.curl(ams)(point), ng.curl(direct)(point), atol=1e-8)
    stats = ssn.AMSCoarseStats()
    assert stats["n_free"] == mesh.nedge - 3 and stats["n_extra"] == 0


def _two_tets(periodic):
    from netgen.csg import Pnt
    from netgen.meshing import (
        Element2D, Element3D, FaceDescriptor, IdentificationType, Mesh, MeshPoint,
    )

    native = Mesh(dim=3)
    points = [
        native.Add(MeshPoint(Pnt(*xyz)))
        for xyz in ((0, 0, 0), (1, 0, 0), (0, 1, 0), (0, 0, 1),
                    (3, 1, 0), (3, 0, 0), (4, 0, 0), (3, 0, 1))
    ]
    for material, indices, name in (
        (1, (0, 1, 2, 3), "air"), (2, (5, 6, 4, 7), "kelvin"),
    ):
        native.SetMaterial(material, name)
        native.Add(Element3D(material, [points[i] for i in indices]))
        boundary = native.Add(
            FaceDescriptor(surfnr=material, domin=material, domout=0, bc=material)
        )
        for face in ((0, 2, 1), (0, 1, 3), (0, 3, 2), (1, 2, 3)):
            native.Add(Element2D(boundary, [points[indices[i]] for i in face]))
    if periodic:
        for master, slave in ((0, 5), (1, 6), (2, 4)):
            native.AddPointIdentification(
                points[master], points[slave], 1, IdentificationType.PERIODIC
            )
    return ng.Mesh(native)


def _pi_constant_rows(mesh, periodic, complex_):
    """Builds the wirebasket AMS and returns Pi * 1 per edge row (x, y, z)."""
    from radia import sparsesolv_ngsolve as ssn

    base = ng.HCurl(mesh, order=2, nograds=True, complex=complex_)
    fes = ng.Periodic(base) if periodic else base
    for dof in range(mesh.nedge, fes.ndof):
        if fes.CouplingType(dof) == ng.COUPLING_TYPE.WIREBASKET_DOF:
            fes.SetCouplingType(dof, ng.COUPLING_TYPE.INTERFACE_DOF)
    u, v = fes.TnT()
    matrix = ng.BilinearForm(fes, symmetric=True)
    matrix += (ng.curl(u) * ng.curl(v) + 1e-6 * u * v) * ng.dx
    pre = ng.Preconditioner(matrix, "bddc", coarsetype="sparsesolv_ams")
    with ng.TaskManager():
        matrix.Assemble()
    assert ssn.AMSCoarseStats()["complex"] == complex_
    ops = ssn.AMSCoarseOperators()
    out = []
    for name in ("Pix", "Piy", "Piz"):
        rows, _, vals, height, _ = ops[name]
        acc = np.zeros(height)
        np.add.at(acc, np.array(rows, dtype=int), np.array(vals))
        out.append(acc[: mesh.nedge])
    del pre
    return np.stack(out, axis=1)


# Fixture identification slave -> master (vertex numbers = insertion order).
_SLAVE_TO_MASTER = {5: 0, 6: 1, 4: 2}


def _expected_pi_one(mesh, fes, periodic):
    """Signed unit-field interpolant per edge row, from the fixture alone.

    Free (representative) rows: oriented from the lower to the higher vertex class
    (slave vertices replaced by their masters), as the periodic element orients it.
    Other rows: the ordinary lower-to-higher physical orientation."""
    vmap = _SLAVE_TO_MASTER if periodic else {}
    free = fes.FreeDofs()
    pts = [np.array(v.point) for v in mesh.vertices]
    rows = []
    for e, edge in enumerate(mesh.edges):
        a, b = (v.nr for v in edge.vertices)
        if free[e]:
            ca, cb = vmap.get(a, a), vmap.get(b, b)
            assert ca != cb
            lo, hi = (a, b) if ca < cb else (b, a)
        else:
            lo, hi = min(a, b), max(a, b)
        rows.append(pts[hi] - pts[lo])
    return np.array(rows), np.array([bool(free[e]) for e in range(mesh.nedge)])


@pytest.mark.parametrize("complex_", [False, True])
@pytest.mark.parametrize("periodic", [True, False])
def test_wirebasket_ams_pi_reproduces_constant_fields(periodic, complex_):
    """Pi * 1 on every edge row is the signed physical edge vector (the edge
    interpolant of the unit fields, oriented like the representative row of G).
    A Pi built from the periodic G alone gives the vector across the period on edges
    touching the paired face; a Pi with the opposite orientation flips the sign."""
    mesh = _two_tets(periodic)
    got = _pi_constant_rows(mesh, periodic, complex_)
    base = ng.HCurl(mesh, order=2, nograds=True, complex=complex_)
    fes = ng.Periodic(base) if periodic else base
    want, active = _expected_pi_one(mesh, fes, periodic)
    assert active.sum() == (mesh.nedge - 3 if periodic else mesh.nedge)
    np.testing.assert_allclose(got[active], want[active], atol=1e-12, rtol=0)
