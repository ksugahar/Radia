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
