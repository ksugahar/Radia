"""A small backward error must not relax the declared true-residual gate."""
import numpy as np
import pytest

pytestmark = pytest.mark.usefixtures("ngsolve_taskmanager")


def _h1_system(convection=0.0):
    from netgen.geom2d import unit_square
    from ngsolve import Mesh, H1, BilinearForm, LinearForm, grad, dx, x, CF

    mesh = Mesh(unit_square.GenerateMesh(maxh=.2))
    fes = H1(mesh, order=2, dirichlet=".*")
    u, v = fes.TnT()
    a = BilinearForm(grad(u) * grad(v) * dx + convection * grad(u)[0] * v * dx).Assemble()
    f = LinearForm((1 + x) * v * dx).Assemble()
    free = np.array([fes.FreeDofs()[i] for i in range(fes.ndof)])
    return fes, a, f, free


def test_large_relative_residual_is_rejected_even_at_small_backward_error():
    from radia._residual_gate import check_true_residual

    fes, a, f, free = _h1_system()
    b = f.vec.FV().NumPy().copy()
    x = np.zeros_like(b)
    x[free] = 1e12                               # huge solution, as in a gauged system
    r = np.zeros_like(b)
    r[free] = 1e-7 * np.linalg.norm(b[free]) / np.sqrt(free.sum())   # ||r||/||b|| = 1e-7
    with pytest.raises(RuntimeError, match="true relative residual"):
        check_true_residual(a.mat, r, x, b, free, "test")
    x[free] = 1.0                                # same residual, ordinary solution size
    with pytest.raises(RuntimeError, match="backward error"):
        check_true_residual(a.mat, r, x, b, free, "test")


def test_nonsymmetric_matrix_through_sparsecholesky_is_rejected():
    from ngsolve import GridFunction
    from radia._residual_gate import check_true_residual

    fes, a, f, free = _h1_system(convection=40.0)
    g = GridFunction(fes)
    g.vec.data = a.mat.Inverse(fes.FreeDofs(), inverse="sparsecholesky") * f.vec
    r = f.vec.CreateVector()
    r.data = f.vec - a.mat * g.vec
    with pytest.raises(RuntimeError, match="true relative residual"):
        check_true_residual(a.mat, r.FV().NumPy(), g.vec.FV().NumPy(),
                               f.vec.FV().NumPy(), free, "test")


def test_nonfinite_solution_is_rejected():
    from radia._residual_gate import check_true_residual

    fes, a, f, free = _h1_system()
    b = f.vec.FV().NumPy().copy()
    x = np.full_like(b, np.inf)
    with pytest.raises(RuntimeError):
        check_true_residual(a.mat, np.zeros_like(b), x, b, free, "test")
