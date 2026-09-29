"""Review finding on the mixed total/reduced Omega matching-trace lane, pinned.

The matching-trace lane recorded the linear residual and returned whatever
CG produced, including a solution that hit ``maxiter``. The residual is now
checked on the returned solution and a miss raises.  (The two findings on
the projected Picard lane were retired with that lane on 2026-09-30.)
"""
from __future__ import annotations

import importlib.util
import math
from pathlib import Path

import pytest

ng = pytest.importorskip("ngsolve")

_HELPERS = None


def helpers():
    """The shared fixtures from the main mixed-Omega test module, by path."""
    global _HELPERS
    if _HELPERS is None:
        path = Path(__file__).with_name("test_kelvin_mixed_omega.py")
        spec = importlib.util.spec_from_file_location("kelvin_mixed_omega_tests", path)
        module = importlib.util.module_from_spec(spec)
        assert spec.loader is not None
        spec.loader.exec_module(module)
        _HELPERS = module
    return _HELPERS


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
