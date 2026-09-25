"""Contracts of radia.p1_newton (first-order Newton, closed-form Jacobian)."""
import importlib.util
import sys
from pathlib import Path

import numpy as np
import pytest

ng = pytest.importorskip("ngsolve")
from radia.p1_newton import lowest_order_gradient, solve_p1_newton  # noqa: E402

MU0 = 4e-7 * np.pi
TABLE = [[0.0, 0.0], [100.0, 0.6], [400.0, 1.2], [2000.0, 1.6], [20000.0, 1.9]]
POINTS = np.asarray([[0.0, 0.0, 0.0], [0.05, 0.0, 0.0], [0.0, 0.0, 0.3], [0.5, 0.5, 0.5]])
CURRENT = 1e7 * ng.CF((-2 * ng.y * (ng.x - .5) * (.7 - ng.x), -(1.2 - 2 * ng.x) * (.01 - ng.y * ng.y), 0))


@pytest.fixture(scope="module")
def box_mesh():
    from netgen.occ import Box, Glue, OCCGeometry, Pnt

    iron = Box(Pnt(-0.2, -0.2, -0.2), Pnt(0.2, 0.2, 0.2))
    iron.mat("iron")
    coil = Box(Pnt(0.5, -0.1, -0.1), Pnt(0.7, 0.1, 0.1))
    coil.mat("coil")
    outer = Box(Pnt(-1, -1, -1), Pnt(1, 1, 1))
    outer.faces.name = "outer"
    outer.mat("air")
    return ng.Mesh(OCCGeometry(Glue([iron, coil, outer - iron - coil])).GenerateMesh(maxh=0.35))


def reduced(mesh, **options):
    settings = dict(source_cf=ng.CF((0.0, 0.0, 0.8)), newton_tolerance=1e-9,
                    field_tolerance=1e-6, observation_points=POINTS)
    settings.update(options)
    return solve_p1_newton(mesh, TABLE, **settings)


def test_lowest_order_gradient_is_create_gradient(box_mesh):
    fes = ng.HCurl(box_mesh, order=1, nograds=True, dirichlet="outer")
    reference, _ = fes.CreateGradient()
    fast = lowest_order_gradient(fes)
    for a, b in zip(reference.CSR(), fast.CSR()):
        np.testing.assert_array_equal(np.asarray(a), np.asarray(b))
    with pytest.raises(ValueError, match="one HCurl dof per edge"):
        lowest_order_gradient(ng.HCurl(box_mesh, order=1))


def test_linear_law_matches_the_independent_linear_solve(box_mesh):
    from radia.vector_potential_solver import VectorPotentialSolver

    mu_r = 500.0
    linear_table = [[0.0, 0.0], [10.0, 10.0 * MU0 * mu_r], [1.0e4, 1.0e4 * MU0 * mu_r]]
    result = solve_p1_newton(box_mesh, linear_table, source_cf=ng.CF((0.0, 0.0, 0.8)),
                             newton_tolerance=1e-10, observation_points=POINTS,
                             linear_solver="direct")
    assert result["stats"]["iterations"] <= 2
    vps = VectorPotentialSolver(box_mesh, mu_r_dict={"iron": mu_r}, order=1)
    vps.set_source_cf(ng.CF((0.0, 0.0, 0.8)))
    with ng.TaskManager():
        vps.solve_linear(dirichlet="outer", solver="direct")
    expected = np.asarray([[float(c) for c in vps.get_B()(box_mesh(*p))] for p in POINTS])
    np.testing.assert_allclose(result["observation_B_T"], expected, rtol=1e-6, atol=1e-8)


def test_linear_solvers_agree_on_the_nonlinear_field(box_mesh):
    beta_zero = reduced(box_mesh)
    iccg = reduced(box_mesh, linear_solver="iccg")
    direct = reduced(box_mesh, linear_solver="direct")
    assert beta_zero["stats"]["beta_zero"] and beta_zero["stats"]["gauge_epsilon"] == 0.0
    np.testing.assert_allclose(iccg["observation_B_T"], beta_zero["observation_B_T"], rtol=1e-6, atol=1e-9)
    # The direct route carries the 1e-6 nu0 mass gauge: a physical perturbation of that order.
    np.testing.assert_allclose(direct["observation_B_T"], beta_zero["observation_B_T"], rtol=1e-4, atol=1e-7)
    for result in (beta_zero, iccg, direct):
        assert result["stats"]["converged"]
        assert all(row["linear_relative_residual"] <= row["linear_tolerance"] * (1 + 1e-9)
                   for row in result["stats"]["history"])


def test_matches_the_validation_lane_on_a_uniform_source(box_mesh):
    # With a uniform source the lane's centroid source equals the element mean,
    # so both implement the same discrete equations.
    lane_path = Path(__file__).resolve().parents[1] / "validation_test" / "c_type_p1_ams_box" / "run_p1_box.py"
    spec = importlib.util.spec_from_file_location("p1_newton_lane_check", lane_path)
    lane = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = lane
    spec.loader.exec_module(lane)
    engine = lane.ReducedAP1Box(box_mesh, ng.CF((0.0, 0.0, 0.8)), linear_solver="ams",
                                cg_tolerance=1e-9, cg_max_iterations=1000, ams_num_smooth=1,
                                source_projection_order=2, gauge_epsilon=0.0, ams_beta_zero=True)
    field, stats, _ = engine.run_newton(lane.SoftIronLaw(TABLE), newton_tolerance=1e-9,
                                        tolerance=1e-6, max_iterations=40, max_halvings=6,
                                        observation=POINTS, inexact_linear=True)
    assert stats["converged"]
    np.testing.assert_allclose(reduced(box_mesh)["observation_B_T"], field, rtol=1e-6, atol=1e-9)


def test_total_a_with_a_divergence_free_current(box_mesh):
    options = dict(current_cf=CURRENT, current_materials="coil", newton_tolerance=1e-9,
                   field_tolerance=1e-6, observation_points=POINTS)
    beta_zero = solve_p1_newton(box_mesh, TABLE, **options)
    direct = solve_p1_newton(box_mesh, TABLE, linear_solver="direct", **options)
    assert np.linalg.norm(beta_zero["observation_B_T"]) > 0
    np.testing.assert_allclose(direct["observation_B_T"], beta_zero["observation_B_T"], rtol=1e-4, atol=1e-7)


def test_incompatible_current_fails_loudly_when_ungauged(box_mesh):
    with pytest.raises(ValueError, match="orthogonal to the discrete gradients"):
        solve_p1_newton(box_mesh, TABLE, current_cf=ng.CF((1e6, 0.0, 0.0)), current_materials="coil")


def test_refuses_to_run_inside_a_taskmanager(box_mesh):
    with ng.TaskManager(), pytest.raises(RuntimeError, match="outside ngsolve.TaskManager"):
        reduced(box_mesh)


def test_non_convergence_and_bad_arguments_raise(box_mesh):
    with pytest.raises(RuntimeError, match="did not converge"):
        reduced(box_mesh, max_iterations=1, newton_tolerance=1e-12)
    with pytest.raises(ValueError, match="exactly one"):
        solve_p1_newton(box_mesh, TABLE)
    with pytest.raises(ValueError, match="gauge_epsilon > 0"):
        reduced(box_mesh, linear_solver="direct", gauge_epsilon=0.0)


def test_vector_potential_solver_entry(box_mesh):
    from radia.vector_potential_solver import VectorPotentialSolver

    vps = VectorPotentialSolver(box_mesh, mu_r_dict={"iron": 1000.0}, order=1)
    vps.set_source_cf(ng.CF((0.0, 0.0, 0.8)))
    vps.solve_nonlinear_newton_p1(TABLE, dirichlet="outer", newton_tolerance=1e-9, field_tolerance=1e-6)
    values = np.asarray([[float(c) for c in vps.get_B()(box_mesh(*p))] for p in POINTS])
    np.testing.assert_allclose(values, reduced(box_mesh)["observation_B_T"], rtol=1e-12, atol=1e-14)
    assert vps._last_nonlinear_stats["converged"]
    with pytest.raises(ValueError, match="order=1"):
        VectorPotentialSolver(box_mesh, mu_r_dict={"iron": 1000.0}, order=2).solve_nonlinear_newton_p1(TABLE)
