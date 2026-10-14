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
    r[free] = 1e-5 * np.linalg.norm(b[free]) / np.sqrt(free.sum())   # ||r||/||b|| = 1e-5
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


@pytest.mark.parametrize("relative", [0.0, 1e-8, 0.999e-6, 1e-6])
def test_shared_linear_gate_accepts_up_to_one_e_minus_six(relative):
    from radia._residual_gate import RELATIVE_LIMIT, check_true_residual

    assert RELATIVE_LIMIT == 1e-6
    # No factorization is involved: this checks the acceptance boundary itself.
    assert check_true_residual(None, np.array([relative]), np.ones(1),
                               np.ones(1), np.array([True]), "boundary") == relative


@pytest.mark.parametrize("relative", [1.001e-6, 1e-5, float("nan"), float("inf")])
def test_shared_linear_gate_rejects_excess_and_nonfinite(relative):
    from radia._residual_gate import check_true_residual

    class Identity:
        def CSR(self):
            return np.ones(1), np.zeros(1, dtype=int), np.array([0, 1])

    with pytest.raises(RuntimeError, match="true relative residual"):
        check_true_residual(Identity(), np.array([relative]), np.ones(1),
                            np.ones(1), np.array([True]), "boundary")


def test_fixed_physical_scale_handles_tiny_correction_rhs_without_backward_admission():
    from radia._residual_gate import check_true_residual

    class Identity:
        def CSR(self):
            return np.ones(1), np.zeros(1, dtype=int), np.array([0, 1])

    free = np.array([True])
    rhs = np.array([1e-20])
    x = np.array([1e12])  # the solution must never enlarge the acceptance scale
    assert check_true_residual(Identity(), np.array([5e-7]), x, rhs, free,
                               "fixed load", reference_norm=1.0) == 5e-7
    with pytest.raises(RuntimeError):
        check_true_residual(Identity(), np.array([2e-6]), x, rhs, free,
                            "fixed load", reference_norm=1.0)
    with pytest.raises(RuntimeError):
        check_true_residual(Identity(), np.array([5e-7]), x, rhs, free, "no scale")


@pytest.mark.parametrize("scale", [-1.0, float("nan"), float("inf")])
def test_physical_reference_scale_must_be_finite_nonnegative(scale):
    from radia._residual_gate import residual_scale
    with pytest.raises(ValueError):
        residual_scale(1.0, scale)
