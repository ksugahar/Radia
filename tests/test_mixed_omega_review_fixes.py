"""Three review findings on the mixed total/reduced Omega lanes, each pinned.

1. The higher-order projected lane capped the permeability at the larger of
   the tabulated node secants and the initial guess. The constitutive law is
   the PCHIP interpolant, whose secant exceeds every node secant between and
   below the nodes: on the Picard test table the nodes top out at 2000 and
   the law reaches about 2620, so the cap rewrote the material by up to 31%,
   and by a different amount for each ``mu_r_initial``. The cap is now taken
   from the law's dense maximum and raised to the law's own target whenever
   the projection asks for more.

2. The same lane declared convergence on the iteration-to-iteration change of
   the projected |B| alone. A stalled update -- small relaxation, or the cap
   above -- makes that change tiny while the material is still far from the
   law; that was reported as ``converged=True``. Convergence now also needs
   the material fixed-point residual, and the final re-solve is checked
   against the law before it is returned.

3. The matching-trace lane recorded the linear residual and returned whatever
   CG produced, including a solution that hit ``maxiter``. The residual is
   now checked on the returned solution and a miss raises.
"""
from __future__ import annotations

import importlib.util
import math
from pathlib import Path

import numpy as np
import pytest

ng = pytest.importorskip("ngsolve")

MU_0 = 4.0e-7 * math.pi
_HELPERS = None


def helpers():
    """The Picard fixtures from the main mixed-Omega test module, by path."""
    global _HELPERS
    if _HELPERS is None:
        path = Path(__file__).with_name("test_kelvin_mixed_omega.py")
        spec = importlib.util.spec_from_file_location("kelvin_mixed_omega_tests", path)
        module = importlib.util.module_from_spec(spec)
        assert spec.loader is not None
        spec.loader.exec_module(module)
        _HELPERS = module
    return _HELPERS


def _law_secant_maximum(bh_table):
    from radia.scalar_potential_solver import _build_bh_interpolator

    table = np.asarray(bh_table, dtype=float)
    positive = table[:, 0] > 0
    law = _build_bh_interpolator(table)
    h = np.geomspace(table[positive, 0].min() * 1.0e-3, table[positive, 0].max(), 20001)
    return float(np.max(np.asarray([law(v) for v in h]) / (MU_0 * h)))


def test_the_projected_cap_comes_from_the_law_not_the_nodes():
    """Finding 1: the ceiling is at least the interpolant's own maximum."""
    h = helpers()
    mesh, h_source, potential, bh_table = h._picard_case()
    table = np.asarray(bh_table, dtype=float)
    node_max = float(np.max(table[1:, 1] / (MU_0 * table[1:, 0])))
    law_max = _law_secant_maximum(bh_table)
    assert law_max > node_max * 1.2, "the table no longer exhibits the defect"

    result = h._picard_solve(
        mesh, h_source, potential, bh_table, order=2, material_update_order=1,
        anderson_depth=0, mu_r_initial=1000.0)
    stats = result["nonlinear_stats"]
    assert stats["converged"]
    lower, upper = stats["physical_permeability_bounds"]
    assert lower == 1.0
    # The solver scans the law on its own grid; this test on a finer one. They
    # agree to ~1e-8, and the bound's job is to stop percent-level clipping,
    # not to reproduce a sampling grid.
    assert upper >= law_max * (1.0 - 1.0e-6), (upper, law_max)
    # The live bound is exp(log(initial)) at minimum; allow the round trip.
    assert upper >= stats["physical_permeability_bounds_initial"][1] * (1.0 - 1.0e-12)
    ceilings = [row["mu_r_upper"] for row in stats["history"]]
    assert ceilings == sorted(ceilings), "the ceiling may only be raised"


def test_the_effective_law_no_longer_depends_on_the_initial_guess():
    """Finding 1, the consequence: two starts, one material, one field."""
    h = helpers()
    mesh, h_source, potential, bh_table = h._picard_case()
    observation = np.array([[0.5, 0.1, 0.2], [0.7, -0.2, 0.1]])
    fields = []
    for initial in (1000.0, 3000.0):
        result = h._picard_solve(
            mesh, h_source, potential, bh_table, order=2, material_update_order=1,
            anderson_depth=0, mu_r_initial=initial, observation_points=observation,
            tolerance=1.0e-8, max_iterations=120)
        assert result["nonlinear_stats"]["converged"]
        fields.append(np.asarray(result["nonlinear_stats"]["observation_field_T"]))
    np.testing.assert_allclose(fields[0], fields[1], rtol=2.0e-6, atol=0.0)


def test_convergence_requires_the_material_residual_too():
    """Finding 2: a converged solve is self-consistent with the law."""
    h = helpers()
    mesh, h_source, potential, bh_table = h._picard_case()
    tolerance = 1.0e-6
    result = h._picard_solve(
        mesh, h_source, potential, bh_table, order=2, material_update_order=1,
        anderson_depth=0, tolerance=tolerance)
    stats = result["nonlinear_stats"]
    assert stats["converged"]
    assert stats["relative_constitutive_change"] <= tolerance
    assert stats["final_state_constitutive_residual"] is not None
    assert stats["final_state_constitutive_residual"] <= tolerance
    assert stats["final_material_state_resolved"]
    assert all("relative_constitutive_change" in row for row in stats["history"])


def test_a_stalled_update_is_not_reported_as_converged():
    """Finding 2, the failure it existed for: tiny relaxation, moving nothing."""
    from radia.kelvin_solver import MixedOmegaPicardNotConverged

    h = helpers()
    mesh, h_source, potential, bh_table = h._picard_case()
    with pytest.raises(MixedOmegaPicardNotConverged, match="did not converge") as info:
        h._picard_solve(
            mesh, h_source, potential, bh_table, order=2, material_update_order=1,
            anderson_depth=0, relaxation=1.0e-3, max_iterations=3, tolerance=1.0e-6)
    stats = info.value.state["nonlinear_stats"]
    assert stats["converged"] is False
    # The material is still far from the law, whatever |B| did.
    assert stats["relative_constitutive_change"] > 1.0e-6
    assert "relative_constitutive_change" in str(info.value)


def _matching_trace_case():
    h = helpers()
    mesh = h._two_region_mesh(maxh=0.8)
    trace_space = ng.H1(mesh, order=2, definedon=mesh.Boundaries("source_total_interface"))
    trace = ng.GridFunction(trace_space)
    trace.vec[:] = 0.0
    source_h = ng.CoefficientFunction((ng.x * ng.x, ng.y, ng.z))
    args = dict(
        mu_r_by_material={"reduced": 1.0, "total": 1.0},
        reduced_materials=("reduced",), total_materials=("total",),
        interface_boundary="source_total_interface", order=2,
        dirichlet_boundary="outer", solver="cg", cg_tolerance=1.0e-10)
    return mesh, source_h, trace, args


def test_a_cg_solve_that_ran_out_of_iterations_raises():
    """Finding 3: hitting maxiter is a failure, not a result."""
    from radia.kelvin_solver import (
        LinearSolveNotConverged,
        solve_magnetostatic_matching_trace_total_reduced_omega,
    )

    mesh, source_h, trace, args = _matching_trace_case()
    with ng.TaskManager():
        with pytest.raises(LinearSolveNotConverged, match="did not converge") as info:
            solve_magnetostatic_matching_trace_total_reduced_omega(
                mesh, source_h, trace, cg_max_iterations=1, **args)
    residual = info.value.residual
    assert residual["free_dofs"]["relative"] > residual["accepted_below"]
    assert "iterations=1" in str(info.value)


def test_a_converged_cg_solve_records_the_bar_it_cleared():
    """Finding 3, the other direction: a real solve still passes, and says how."""
    from radia.kelvin_solver import solve_magnetostatic_matching_trace_total_reduced_omega

    mesh, source_h, trace, args = _matching_trace_case()
    with ng.TaskManager():
        result = solve_magnetostatic_matching_trace_total_reduced_omega(
            mesh, source_h, trace, **args)
    residual = result["linear_residual"]
    assert math.isfinite(result["linear_residual_relative"])
    assert result["linear_residual_relative"] <= residual["accepted_below"]
    assert residual["accepted_below"] == pytest.approx(100.0 * 1.0e-10)
