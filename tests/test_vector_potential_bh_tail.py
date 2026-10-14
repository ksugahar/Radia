"""High-field energy and field consistency for reduced-A Newton solves."""

import numpy as np
import pytest

ng = pytest.importorskip("ngsolve")
from ngsolve.meshes import MakeStructured3DMesh

from radia.vector_potential_solver import MU_0, VectorPotentialSolver, _coercive_bh_law


def test_bh_tail_preserves_energy_derivative_and_endpoint_continuity():
    h, energy, tail = _coercive_bh_law(
        [0.0, 0.2, 1.0, 3.0], [0.0, 100.0, 1000.0, 1.0e6])
    with ng.TaskManager():
        mesh = MakeStructured3DMesh(hexes=False, nx=1, ny=1, nz=1)
        volume = ng.Integrate(ng.CF(1.0), mesh)

        def evaluate(law, b):
            return float(ng.Integrate(law(ng.CF(b)), mesh) / volume)

        endpoint = tail["B_max_ext"]
        delta = 1.0e-5
        for b in (endpoint, endpoint + 1.0, endpoint + 20.0):
            derivative = (evaluate(energy, b + delta)
                          - evaluate(energy, b - delta)) / (2 * delta)
            assert derivative == pytest.approx(evaluate(h, b), rel=2e-8)
        assert evaluate(h, endpoint + delta) - evaluate(h, endpoint - delta) == (
            pytest.approx(2 * delta / MU_0, rel=2e-8))
        assert evaluate(energy, endpoint + 20.0) > evaluate(energy, endpoint) > 0


@pytest.mark.parametrize("source_b", [10.0, 20.0])
def test_newton_postprocessing_keeps_high_field_constitutive_response(source_b):
    with ng.TaskManager():
        mesh = MakeStructured3DMesh(hexes=False, nx=1, ny=1, nz=1)
        material = mesh.GetMaterials()[0]
        solver = VectorPotentialSolver(mesh, iron_domains=material, order=1)
        solver.set_source_cf(ng.CF((source_b, 0.0, 0.0)))
        solver.solve_nonlinear_newton(
            [[0.0, 0.0], [100.0, 0.2], [1000.0, 1.0], [1.0e6, 3.0]],
            solver="direct", maxiter=3, tol=1e-10, verbose=False)
        point = mesh(0.23, 0.31, 0.41)
        observed_b = np.asarray(solver.get_B()(point))
        observed_h = np.asarray(solver.get_H()(point))
    np.testing.assert_allclose(observed_b, [source_b, 0.0, 0.0], atol=1e-9)
    expected_h = 1.0e6 + (source_b - 3.0) / MU_0
    np.testing.assert_allclose(observed_h, [expected_h, 0.0, 0.0], rtol=1e-10, atol=1e-7)
