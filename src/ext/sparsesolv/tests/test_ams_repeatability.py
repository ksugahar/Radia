"""Nearly singular gradient correction must not depend on task scheduling."""
import numpy as np
import pytest
import ngsolve as ng
from netgen.csg import unit_cube
from radia.sparsesolv_ngsolve import HypreBasedAMSPreconditioner


@pytest.mark.parametrize("beta", [1.0, 1e-6])
@pytest.mark.parametrize("subspace_solver", [0, 1])
def test_parallel_gradient_correction_repeatability(beta, subspace_solver):
    ng.SetNumThreads(1)
    mesh = ng.Mesh(unit_cube.GenerateMesh(maxh=0.35))
    fes = ng.HCurl(mesh, order=1, nograds=True, dirichlet=".*")
    u, v = fes.TnT()
    a = ng.BilinearForm(fes)
    a += (ng.IfPos(ng.x-0.5, 1000, 1)*ng.curl(u)*ng.curl(v) + beta*u*v)*ng.dx
    a.Assemble()
    grad, _ = fes.CreateGradient()
    coords = np.asarray(mesh.ngmesh.Coordinates())
    pre = HypreBasedAMSPreconditioner(a.mat, grad, freedofs=fes.FreeDofs(),
        coord_x=coords[:,0].tolist(), coord_y=coords[:,1].tolist(),
        coord_z=coords[:,2].tolist(), subspace_solver=subspace_solver)
    exact, rhs, out = a.mat.CreateColVector(), a.mat.CreateRowVector(), a.mat.CreateColVector()
    exact.FV().NumPy()[:] = np.random.default_rng(1004).normal(size=fes.ndof)*np.asarray(list(fes.FreeDofs()))
    rhs.data = a.mat*exact
    reference = None
    try:
        for threads in (1, 4):
            ng.SetNumThreads(threads)
            with ng.TaskManager():
                for _ in range(12):
                    pre.Mult(rhs, out)
                    value = out.FV().NumPy().copy()
                    if reference is None:
                        reference = value
                    else:
                        np.testing.assert_array_equal(value, reference)
    finally:
        ng.SetNumThreads(1)
