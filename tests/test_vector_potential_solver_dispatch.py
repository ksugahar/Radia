import importlib.util
import inspect
import sys
import types
from pathlib import Path

import numpy as np
import pytest


MODULE_PATH = (
    Path(__file__).resolve().parents[1]
    / "src"
    / "radia"
    / "vector_potential_solver.py"
)
KELVIN_SOLVER_PATH = MODULE_PATH.with_name("kelvin_solver.py")
SPEC = importlib.util.spec_from_file_location("vector_potential_solver_under_test", MODULE_PATH)
MODULE = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
SPEC.loader.exec_module(MODULE)
VectorPotentialSolver = MODULE.VectorPotentialSolver


def selector(order):
    value = VectorPotentialSolver.__new__(VectorPotentialSolver)
    value.order = order
    value._kelvin_region = None
    return value


def test_high_order_auto_dispatches_to_bddc_before_ams_setup():
    assert selector(2)._select_solver(200_001, "auto") == "bddc"


def test_high_order_explicit_ams_fails_loud():
    with pytest.raises(ValueError, match="requires HCurl order=1"):
        selector(2)._select_solver(200_001, "ams")


def test_order_one_auto_uses_current_ams_factory(monkeypatch):
    package = types.ModuleType("radia")
    package.__path__ = []
    ssn = types.ModuleType("radia.sparsesolv_ngsolve")
    ssn.HypreBasedAMSPreconditioner = object()
    package.sparsesolv_ngsolve = ssn
    monkeypatch.setitem(sys.modules, "radia", package)
    monkeypatch.setitem(sys.modules, "radia.sparsesolv_ngsolve", ssn)
    assert selector(1)._select_solver(200_001, "auto") == "ams"


def test_periodic_kelvin_auto_uses_bddc_for_large_system():
    value = selector(1)
    value._kelvin_region = "kelvin"
    assert value._select_solver(1_000_000, "auto") == "bddc"


def test_periodic_kelvin_auto_keeps_direct_for_small_system():
    value = selector(1)
    value._kelvin_region = "kelvin"
    assert value._select_solver(200_000, "auto") == "direct"


def test_periodic_kelvin_accepts_bddc():
    value = selector(2)
    value._kelvin_region = "kelvin"
    assert value._select_solver(1_000_000, "bddc") == "bddc"


def test_periodic_kelvin_rejects_ams_low_order_auxiliary_space():
    value = selector(1)
    value._kelvin_region = "kelvin"
    with pytest.raises(ValueError, match="Periodic low-order coupling"):
        value._select_solver(50_000, "ams")


def test_gauge_defaults_are_stable_and_shared_across_solve_paths():
    assert MODULE.DEFAULT_GAUGE_EPSILON == pytest.approx(1.0e-6)
    for method_name in (
        "solve_linear", "solve_nonlinear", "solve_nonlinear_newton"
    ):
        signature = inspect.signature(getattr(VectorPotentialSolver, method_name))
        assert signature.parameters["eps"].default == pytest.approx(1.0e-6)


def test_kelvin_gauge_defaults_to_physical_value_but_can_be_split():
    value = selector(2)
    value._kelvin_region = "kelvin"
    assert value._resolve_gauge_epsilons(1.0e-6, None) == (
        pytest.approx(1.0e-6), pytest.approx(1.0e-6))
    assert value._resolve_gauge_epsilons(2.0e-6, 5.0e-6) == (
        pytest.approx(2.0e-6), pytest.approx(5.0e-6))


@pytest.mark.parametrize("eps", [0.0, -1.0, float("inf"), float("nan")])
def test_gauge_rejects_non_positive_or_non_finite_values(eps):
    value = selector(2)
    with pytest.raises(ValueError, match="positive finite"):
        value._resolve_gauge_epsilons(eps, None)


def test_nonlinear_reduced_a_assembles_split_gauge_coefficients():
    source = MODULE_PATH.read_text(encoding="utf-8")
    nonlinear = source.split("def solve_nonlinear(", 1)[1].split(
        "def solve_hysteresis(", 1
    )[0]

    assert "physical_gauge_coeff = physical_eps * nu_air" in nonlinear
    assert "kelvin_gauge_coeff" in nonlinear
    assert "physical_gauge_coeff * InnerProduct" in nonlinear
    assert "kelvin_gauge_coeff * InnerProduct" in nonlinear


def test_newton_reduced_a_uses_vacuum_scaled_gauge_energy():
    source = MODULE_PATH.read_text(encoding="utf-8")
    newton = source.split("def solve_nonlinear_newton(", 1)[1].split(
        "def solve_nonlinear(", 1
    )[0]

    assert "eps, _ = self._resolve_gauge_epsilons(eps, None)" in newton
    assert "gauge_coeff = eps * nu_0" in newton
    assert "gauge_coeff / 2.0 * InnerProduct(A, A)" in newton


def test_linear_solve_has_finite_solution_and_true_residual_gate():
    source = MODULE_PATH.read_text(encoding="utf-8")
    linear = source.split("def solve_linear(", 1)[1].split(
        "def solve_nonlinear_newton(", 1
    )[0]

    assert "np.all(np.isfinite(solution_values))" in linear
    assert "linear_residual.data = f.vec - a.mat * self._A_gf.vec" in linear
    assert "_relative_residual_on_free_dofs" in linear
    assert "LINEAR_RELATIVE_RESIDUAL_LIMIT" in linear


def test_relative_residual_ignores_dirichlet_reaction_and_rhs():
    matrix = np.diag([2.0, 3.0, 4.0])
    solution = np.array([0.5, 2.0 / 3.0, 0.0])
    rhs = np.array([1.0, 2.0, 7.0])
    residual = rhs - matrix @ solution
    free = np.array([True, True, False])

    assert np.linalg.norm(residual) > 0.0
    assert MODULE._relative_residual_on_free_dofs(
        residual, rhs, free) == pytest.approx(0.0)


def test_low_level_kelvin_a_helpers_use_the_stable_gauge_default():
    source = KELVIN_SOLVER_PATH.read_text(encoding="utf-8")
    assert "gauge_eps=1e-8" not in source
    assert source.count("gauge_eps=1e-6") == 3
