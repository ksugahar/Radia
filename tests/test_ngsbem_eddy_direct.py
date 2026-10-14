"""Production scalar eddy solves retain their operator with SparseCholesky."""
import numpy as np
import pytest

ngs = pytest.importorskip("ngsolve")

from radia.ngsbem_eddy import EddyCurrentFEMBEM, create_conductor_mesh


@pytest.mark.parametrize("frequency", [0.0, 5000.0])
@pytest.mark.parametrize("order", [2, 3])
def test_dirichlet_solve_matches_independent_dense_system(frequency, order):
    from scipy.linalg import solve

    mesh = create_conductor_mesh(0.01, 0.01, 0.004, maxh=0.005)
    solver = EddyCurrentFEMBEM(mesh, order=order)
    with ngs.TaskManager():
        solver.assemble_fem(frequency)
        solver.set_source_field(1 + 20 * ngs.x + 10j * ngs.y)
        solution = solver.solve()
    # Assemble the independent free-row Dirichlet problem from native COO,
    # retaining the returned boundary values and all complex coefficients.
    from scipy.sparse import coo_matrix
    rows, cols, data = solver._a.mat.COO()
    matrix = coo_matrix((np.asarray(data), (rows, cols)),
                        shape=(solver._fes.ndof, solver._fes.ndof)).toarray()
    free = np.array(list(solver._fes.FreeDofs()), dtype=bool)
    values = solution.vec.FV().NumPy().copy()
    rhs = -matrix[np.ix_(free, ~free)] @ values[~free]
    reference = solve(matrix[np.ix_(free, free)], rhs)
    np.testing.assert_allclose(values[free], reference, rtol=1e-10, atol=1e-11)
    assert solver.true_relative_residual < 1e-10


@pytest.mark.parametrize("frequency", [500.0, 5000.0])
@pytest.mark.parametrize("amplitude", [1e-3, 1e-9])
@pytest.mark.parametrize("order", [1, 2])
def test_fembem_sparsecholesky_blocks_converge(frequency, amplitude, order):
    mesh = create_conductor_mesh(0.01, 0.01, 0.004, maxh=0.006)
    solver = EddyCurrentFEMBEM(mesh, order=order)
    with ngs.TaskManager():
        solver.assemble_fembem(frequency, intorder=6)
        solution = solver.solve(B_ext=[0, 0, amplitude], mode="fembem")
    assert np.all(np.isfinite(solution.vec.FV().NumPy()))
    assert solver.true_relative_residual <= 1e-8


def test_fembem_rejects_unconverged_iterate(monkeypatch):
    import ngsolve.solvers

    def zero_iterate(*, b, **kwargs):
        result = b.CreateVector()
        result[:] = 0
        return result

    monkeypatch.setattr(ngsolve.solvers, "GMRes", zero_iterate)
    mesh = create_conductor_mesh(0.01, 0.01, 0.004, maxh=0.006)
    solver = EddyCurrentFEMBEM(mesh, order=1)
    with ngs.TaskManager():
        solver.assemble_fembem(500.0, intorder=6)
        with pytest.raises(RuntimeError, match="FEM-BEM true relative residual"):
            solver.solve(B_ext=[0, 0, 1e-3], mode="fembem")
    assert not solver._solved
