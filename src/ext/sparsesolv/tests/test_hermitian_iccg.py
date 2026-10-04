"""Production Hermitian ICCG: independent dense oracle and rejection checks."""
import numpy as np
import pytest
import ngsolve as ng
from radia.sparsesolv_ngsolve import SparseSolvSolver


def system(a, values):
    n = len(a)
    from netgen.geom2d import unit_square
    mesh = ng.Mesh(unit_square.GenerateMesh(maxh=1))
    space = ng.FESpace([ng.NumberSpace(mesh, complex=True) for _ in range(n)], complex=True)
    u, v = space.TnT()
    form = ng.BilinearForm(space)
    form += sum(u[i]*v[j] for i in range(n) for j in range(n))*ng.dx
    form.Assemble()
    matrix = form.mat
    i, j, _ = matrix.COO()
    matrix.AsVector().FV().NumPy()[:] = a[np.asarray(i), np.asarray(j)]
    rhs, solution = matrix.CreateRowVector(), matrix.CreateColVector()
    rhs.FV().NumPy()[:] = values
    solution.FV().NumPy()[:] = 0
    return matrix, rhs, solution


@pytest.mark.parametrize("threads", [1, 4])
@pytest.mark.parametrize("scaling", [False, True])
@pytest.mark.parametrize("ordering", ["natural", "abmc", "rcm_spmv"])
@pytest.mark.parametrize("hermitian", [False, True])
def test_complex_dense_oracle(threads, scaling, ordering, hermitian):
    rng = np.random.default_rng(42)
    q = rng.normal(size=(19, 19)) + 1j*rng.normal(size=(19, 19))
    a = (q @ q.conj().T if hermitian else q.real @ q.real.T + 1j*(q.imag @ q.imag.T)) + 3*np.eye(19)
    b = rng.normal(size=19) + 1j*rng.normal(size=19)
    matrix, rhs, solution = system(a, b)
    solver = SparseSolvSolver(matrix, conjugate=hermitian, diagonal_scaling=scaling,
        use_abmc=ordering != "natural", abmc_use_rcm=ordering == "rcm_spmv",
        abmc_reorder_spmv=ordering == "rcm_spmv", tol=1e-12, maxiter=100)
    ng.SetNumThreads(threads)
    try:
        with ng.TaskManager():
            result = solver.Solve(rhs, solution)
        assert result.converged
        assert result.iterations == 1  # Full-pattern IC is exact here.
        np.testing.assert_allclose(solution.FV().NumPy(), np.linalg.solve(a, b),
                                   rtol=1e-10, atol=1e-11)
        assert np.linalg.norm(a@solution.FV().NumPy()-b)/np.linalg.norm(b) < 1e-11
    finally:
        ng.SetNumThreads(1)


@pytest.mark.parametrize("a", [
    [[2, 1j], [1j, 2]],  # symmetric but not Hermitian
    [[2+1j, 0], [0, 2]],
    [[complex(1.7e308, 1.7e308), 0], [0, 2]],
    [[2, 1j], [0, 2]],
    [[-2, 0], [0, 2]],
    [[2, complex(float('nan'), 0)], [0, 2]],
])
def test_invalid_hermitian_rejected_even_for_zero_rhs(a):
    matrix, rhs, solution = system(np.asarray(a, complex), [0, 0])
    with pytest.raises(ValueError):
        SparseSolvSolver(matrix, conjugate=True).Solve(rhs, solution)


def test_indefinite_negative_curvature_rejected():
    matrix, rhs, solution = system(np.array([[1, 2], [2, 1]], complex), [1, -1])
    with pytest.raises(RuntimeError, match="curvature"):
        SparseSolvSolver(matrix, conjugate=True, shift=3, auto_shift=False).Solve(rhs, solution)


def test_reused_ic_configuration_matches_fresh_factor():
    from radia.sparsesolv_ngsolve import ICPreconditioner
    rng = np.random.default_rng(5)
    q = rng.normal(size=(13, 13))
    matrix, rhs, output = system((q@q.T + 2*np.eye(13)).astype(complex), np.ones(13))
    reused = ICPreconditioner(matrix)
    for abmc, scaling, block in [(False, True, 4), (True, False, 3),
                                  (True, True, 5), (False, False, 4)]:
        fresh = ICPreconditioner(matrix)
        for pre in (reused, fresh):
            pre.use_abmc = abmc
            pre.diagonal_scaling = scaling
            pre.abmc_block_size = block
            pre.Update()
        output.data = reused*rhs
        expected = rhs.CreateVector(); expected.data = fresh*rhs
        np.testing.assert_array_equal(output.FV().NumPy(), expected.FV().NumPy())
