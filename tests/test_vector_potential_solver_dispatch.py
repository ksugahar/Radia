import importlib.util
import sys
import types
from pathlib import Path

import pytest


MODULE_PATH = (
    Path(__file__).resolve().parents[1]
    / "src"
    / "radia"
    / "vector_potential_solver.py"
)
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


def test_periodic_kelvin_auto_uses_direct_solver():
    value = selector(1)
    value._kelvin_region = "kelvin"
    assert value._select_solver(1_000_000, "auto") == "direct"


@pytest.mark.parametrize("solver", ["ams", "bddc"])
def test_periodic_kelvin_rejects_iterative_solvers_that_produce_nan(solver):
    value = selector(1)
    value._kelvin_region = "kelvin"
    with pytest.raises(ValueError, match="Periodic Kelvin HCurl"):
        value._select_solver(50_000, solver)


def test_nonlinear_reduced_a_defines_dimensionless_gauge_coefficient():
    source = MODULE_PATH.read_text(encoding="utf-8")
    nonlinear = source.split("def solve_nonlinear(", 1)[1].split(
        "def solve_hysteresis(", 1
    )[0]

    assert "gauge_coeff = eps * nu_air" in nonlinear
    assert "gauge_coeff * InnerProduct(A_trial, v)" in nonlinear


def test_newton_reduced_a_uses_vacuum_scaled_gauge_energy():
    source = MODULE_PATH.read_text(encoding="utf-8")
    newton = source.split("def solve_nonlinear_newton(", 1)[1].split(
        "def solve_nonlinear(", 1
    )[0]

    assert "gauge_coeff = eps * nu_0" in newton
    assert "gauge_coeff / 2.0 * InnerProduct(A, A)" in newton
