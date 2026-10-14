"""AMS preconditions the Periodic (Kelvin exterior) order-1 reduced-A system.

``CreateGradient`` has no Periodic variant, so the solver builds the discrete
gradient of the Periodic H1 space itself.  The gradient must be exact, every row
an edge-vertex pair, and AMS-CG must reproduce the direct solution.  Measured on
validation runtime (NGSolve 6.2.2606): 23/27/31 CG iterations for 4k/12k/48k tetrahedra, curl A
within 3e-10 of SparseCholesky.
"""
import math

import numpy as np
import pytest

ng = pytest.importorskip("ngsolve")
pytest.importorskip("radia.sparsesolv_ngsolve")

OFFSET = (3.0, 0.0, 0.0)


@pytest.fixture(scope="module")
def kelvin_mesh():
    from netgen.meshing import IdentificationType
    from netgen.occ import Glue, OCCGeometry, Pnt, Sphere

    iron = Sphere(Pnt(0, 0, 0), 0.35)
    iron.mat("iron")
    iron.maxh = 0.10
    shell = Sphere(Pnt(0, 0, 0), 1.0)
    shell.mat("air")
    for face in shell.faces:
        face.name = "kelvin_int"
    air = shell - iron
    air.mat("air")
    exterior = Sphere(Pnt(*OFFSET), 1.0)
    exterior.mat("kelvin")
    for face in exterior.faces:
        face.name = "kelvin_ext"
    for vertex in exterior.vertices:
        if vertex.p[2] > 0.0:
            vertex.name = "GND"
    geometry = Glue([air, iron, exterior])
    inner = outer = None
    for solid in geometry.solids:
        for face in solid.faces:
            if face.name == "kelvin_int":
                inner = face
            elif face.name == "kelvin_ext":
                outer = face
    inner.Identify(outer, "periodic", IdentificationType.PERIODIC)
    with ng.TaskManager():
        return ng.Mesh(OCCGeometry(geometry).GenerateMesh(maxh=0.26, grading=0.6))


def _solver(mesh):
    from radia.vector_potential_solver import VectorPotentialSolver

    return VectorPotentialSolver(mesh, iron_domains="iron", mu_r=1000.0, order=1,
                                 kelvin_region="kelvin", kelvin_radius=1.0,
                                 kelvin_center=OFFSET)


def test_periodic_gradient_is_exact_and_edge_vertex_paired(kelvin_mesh):
    solver = _solver(kelvin_mesh)
    fes = solver._make_hcurl_space("GND")
    gradient, h1_ndof = solver._periodic_lowest_order_gradient(fes)
    assert h1_ndof == kelvin_mesh.nv
    rows, _, vals = (np.asarray(v) for v in gradient.COO())
    nonzero = vals != 0.0
    assert set(np.abs(vals[nonzero]).tolist()) == {1.0}
    assert np.all(np.bincount(rows[nonzero], minlength=fes.ndof) == 2)
    h1 = ng.Periodic(ng.H1(kelvin_mesh, order=1))
    phi = ng.GridFunction(h1)
    phi.vec.FV().NumPy()[:] = np.random.default_rng(3).standard_normal(h1.ndof)
    field = ng.GridFunction(fes)
    field.vec.data = gradient * phi.vec
    with ng.TaskManager():
        error = ng.Integrate(ng.InnerProduct(field - ng.grad(phi), field - ng.grad(phi)), kelvin_mesh)
        scale = ng.Integrate(ng.InnerProduct(ng.grad(phi), ng.grad(phi)), kelvin_mesh)
    assert math.sqrt(abs(error)/scale) < 1e-12


def test_periodic_kelvin_ams_matches_direct(kelvin_mesh):
    from ngsolve import InnerProduct, curl, dx
    from ngsolve.krylovspace import CGSolver
    from radia.kelvin_source import kelvin_nu_factor_3d_cf

    solver = _solver(kelvin_mesh)
    fes = solver._make_hcurl_space("GND")
    nu_air = 1.0/(4e-7*math.pi)
    nu = kelvin_mesh.MaterialCF({"iron": nu_air/1000.0}, default=nu_air)
    eps = 1e-6*nu_air
    u, v = fes.TnT()
    a = ng.BilinearForm(fes, symmetric=True)
    for material in ("air", "iron"):
        a += nu*InnerProduct(curl(u), curl(v))*dx(material) + eps*InnerProduct(u, v)*dx(material)
    a += (nu_air*kelvin_nu_factor_3d_cf(OFFSET, 1.0)*InnerProduct(curl(u), curl(v))
          + eps*InnerProduct(u, v))*dx("kelvin")
    f = ng.LinearForm(fes)
    f += InnerProduct(ng.CF((0, 0, 1.0))*kelvin_mesh.MaterialCF({"iron": 1.0}, default=0.0),
                      curl(v))*dx
    with ng.TaskManager():
        a.Assemble()
        f.Assemble()
        direct = ng.GridFunction(fes)
        direct.vec.data = a.mat.Inverse(fes.FreeDofs(), inverse="sparsecholesky")*f.vec
    # Native AMS setup must run outside TaskManager; the solve may run inside.
    pre = solver._setup_ams_preconditioner(a.mat, fes, eps)
    with ng.TaskManager():
        cg = CGSolver(a.mat, pre, maxiter=500, tol=1e-10, printrates=False)
        iterative = ng.GridFunction(fes)
        iterative.vec.data = cg*f.vec
        difference = curl(iterative) - curl(direct)
        region = kelvin_mesh.Materials("air|iron")
        error = ng.Integrate(InnerProduct(difference, difference), kelvin_mesh, definedon=region)
        scale = ng.Integrate(InnerProduct(curl(direct), curl(direct)), kelvin_mesh, definedon=region)
    assert solver._select_solver(fes.ndof, "ams") == "ams"
    assert cg.iterations < 60
    assert math.sqrt(abs(error)/scale) < 1e-8
