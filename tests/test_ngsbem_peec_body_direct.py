"""Verify the production multi-source body solve against a dense oracle."""
import numpy as np
import pytest

ngs = pytest.importorskip("ngsolve")


@pytest.mark.parametrize("frequency", [500.0, 5000.0])
@pytest.mark.parametrize("order", [2, 3])
def test_body_coupling_preserves_free_equations(monkeypatch, frequency, order):
    from scipy.linalg import solve
    from scipy.sparse import coo_matrix
    import radia.ngsbem_peec_body_coupling as body
    from radia.ngsbem_eddy import create_conductor_mesh

    actual = body._apply_dirichlet_inverse
    residuals = []

    def checked(matrix, inverse, fes, field):
        rows, cols, data = matrix.COO()
        dense = coo_matrix((np.asarray(data), (rows, cols)),
                           shape=(fes.ndof, fes.ndof)).toarray()
        free = np.array(list(fes.FreeDofs()), dtype=bool)
        boundary = field.vec.FV().NumPy().copy()
        rhs = -dense[np.ix_(free, ~free)] @ boundary[~free]
        expected = solve(dense[np.ix_(free, free)], rhs)
        residual = actual(matrix, inverse, fes, field)
        np.testing.assert_allclose(field.vec.FV().NumPy()[free], expected,
                                   rtol=1e-10, atol=1e-11)
        np.testing.assert_array_equal(field.vec.FV().NumPy()[~free], boundary[~free])
        residuals.append(residual)
        return residual

    monkeypatch.setattr(body, "_apply_dirichlet_inverse", checked)
    mesh = create_conductor_mesh(0.01, 0.01, 0.004, maxh=0.005)
    coupled = body.CoupledPEECBody(mesh, body_sigma=5.8e7, fem_order=order)
    for z in (0.01, -0.01):
        coupled.add_filament_group(body.create_rectangular_loop_segments(
            (0, 0, z), 0.03, 0.03))
    with ngs.TaskManager():
        impedance = coupled.compute_delta_Z(frequency)
        permeability = coupled.compute_mu_eff(frequency)
        dc = coupled.compute_delta_Z(0)
    assert impedance.shape == (2, 2)
    assert np.all(np.isfinite(impedance))
    assert np.isfinite(permeability)
    assert len(residuals) == 3
    assert max(residuals) < 1e-10
    np.testing.assert_array_equal(dc, 0)
