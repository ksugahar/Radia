"""Reduced-A Picard loop: history, per-element warm start, constrained Anderson mixing."""
import math

import numpy as np
import pytest

ng = pytest.importorskip("ngsolve")
pytest.importorskip("netgen.occ")

from radia.vector_potential_solver import VectorPotentialSolver  # noqa: E402

MU_0 = 4.0e-7 * math.pi


def _iron_in_air(maxh=0.45):
    from netgen.occ import Box, Glue, OCCGeometry, Pnt

    iron = Box(Pnt(-0.4, -0.4, -0.4), Pnt(0.4, 0.4, 0.4))
    iron.mat("iron")
    outer = Box(Pnt(-1.5, -1.5, -1.5), Pnt(1.5, 1.5, 1.5))
    # Name the outer boundary before the boolean so the iron/air interface
    # keeps its own label: a Dirichlet interface would isolate the iron and turn
    # the Picard loop into a two-iteration no-op.
    for face in outer.faces:
        face.name = "outer"
    air = outer - iron
    air.mat("air")
    with ng.TaskManager():
        mesh = ng.Mesh(OCCGeometry(Glue([iron, air])).GenerateMesh(maxh=maxh))
    boundaries = mesh.GetBoundaries()
    assert boundaries.count("outer") == 6, boundaries
    assert len(boundaries) > 6, boundaries
    return mesh


def _bh_table():
    return [[0.0, 0.0], [100.0, 0.6], [400.0, 1.2], [2000.0, 1.6], [20000.0, 1.9]]


def _solve(mesh, **overrides):
    solver = VectorPotentialSolver(mesh, iron_domains="iron", mu_r=1000.0, order=1)
    # A source that drives the block into the knee of the table, with a tolerance
    # tight enough that the damped Picard loop needs a genuine iteration count.
    solver.set_source_cf(ng.CoefficientFunction((0.0, 0.0, 1.3)))
    settings = dict(tol=1.0e-6, maxiter=80, relax=0.3, dirichlet="outer",
                    verbose=False, solver="direct")
    settings.update(overrides)
    with ng.TaskManager():
        solver.solve_nonlinear(_bh_table(), **settings)
    return solver, dict(solver._last_nonlinear_stats)


def test_reduced_a_picard_records_history_and_per_element_state():
    mesh = _iron_in_air()
    observation = np.array([[0.0, 0.0, 0.9], [0.5, 0.2, 0.0]])
    _, stats = _solve(mesh, observation_points=observation)
    assert stats["converged"]
    assert stats["anderson_depth"] == 0
    assert stats["warm_start"] is False
    assert stats["relaxation"] == 0.3
    assert len(stats["history"]) == stats["iterations"]
    assert "relative_B_change" not in stats["history"][0]
    assert all("relative_B_change" in row for row in stats["history"][1:])
    assert stats["history"][-1]["relative_B_change"] == stats["final_relative_change"]
    assert all("observation_relative_change" in row for row in stats["history"][1:])
    assert np.asarray(stats["observation_field_T"]).shape == (2, 3)
    iron = [el.nr for el in mesh.Elements(ng.VOL) if str(el.mat) == "iron"]
    assert stats["element_numbers"] == iron
    nu = np.asarray(stats["nu_elements"])
    assert nu.shape == (len(iron),)
    assert np.all(nu > 0.0) and np.all(nu <= 1.0 / MU_0)


def test_reduced_a_picard_warm_start_and_anderson_reach_the_same_field():
    mesh = _iron_in_air()
    cold_solver, cold = _solve(mesh)
    warm_solver, warm = _solve(mesh, nu_initial=np.asarray(cold["nu_elements"]))
    assert warm["warm_start"] is True
    assert warm["converged"]
    assert cold["iterations"] >= 5
    assert warm["iterations"] <= 4
    assert warm["iterations"] < cold["iterations"]
    mixed_solver, mixed = _solve(mesh, anderson_depth=2)
    assert mixed["converged"]
    assert mixed["anderson"]["accelerated_steps"] >= 1
    assert mixed["iterations"] <= cold["iterations"]
    probe = mesh(0.0, 0.0, 0.9)
    reference = np.asarray(cold_solver.get_B()(probe))
    np.testing.assert_allclose(np.asarray(warm_solver.get_B()(probe)), reference, rtol=2.0e-3)
    np.testing.assert_allclose(np.asarray(mixed_solver.get_B()(probe)), reference, rtol=5.0e-3)
    with pytest.raises(ValueError, match="one reluctivity per iron element"):
        _solve(mesh, nu_initial=np.ones(2))
