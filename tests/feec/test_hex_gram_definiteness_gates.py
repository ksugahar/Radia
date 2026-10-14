"""The HEX Gram definiteness gates must not pass an unconverged computation.

Two false-pass routes existed and are locked shut here.

1. ``floor_scan`` treated "the C++ returned without raising" as convergence.
   The production CG returns normally when it runs out of iterations; only a
   breakdown raises.  The convergence flag and the true relative residual are
   already reported in the solve timings, so the gate must read them.

2. The spectral gate accepted whatever LOBPCG returned.  A Rayleigh-Ritz value
   bounds the smallest eigenvalue from ABOVE, so an unconverged smallest Ritz
   value clearing the floor is not evidence that no negative eigenvalue exists.
   Non-convergence must withhold judgement.

These run on test doubles: the point is the decision logic, not the kernel.
"""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import numpy as np
import pytest

REPO = next(p for p in Path(__file__).resolve().parents if (p / "src" / "radia").exists())
VALIDATOR = (REPO / "validation_test" / "esrf_three_engine"
             / "validate_hex_gram_definiteness.py")


def _load_validator():
    if not VALIDATOR.is_file():
        pytest.skip(f"validator not present: {VALIDATOR}")
    spec = importlib.util.spec_from_file_location("_hex_gram_validator", VALIDATOR)
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


class _Solver:
    """Stands in for RadHACApKChargeGram's linear-material solve."""

    def __init__(self, converged, residual, iters=1200):
        self._converged = converged
        self._residual = residual
        self._iters = iters
        self.timings = True

    def solve_configured_linear_material_auto_prec(self, floor, rhs, tol, maxit):
        result = {"iters": self._iters, "prec_min": 1.0, "prec_max": 2.0}
        if self.timings:
            result["timings"] = {
                "last_solve_converged": 1.0 if self._converged else 0.0,
                "last_solve_final_relative_residual": self._residual,
            }
        return result


def test_floor_scan_rejects_a_solve_that_ran_out_of_iterations():
    validator = _load_validator()
    solver = _Solver(converged=False, residual=4.2e-3)
    scan = validator.floor_scan(solver, np.zeros(4), 1.0e3, [1.0, 2.0], 1.0e-8, 1200)
    assert scan[0]["converged"] is False
    assert scan[0]["relative_residual"] == pytest.approx(4.2e-3)
    # The scan stops at the first failure, so the caller's length check fires too.
    assert len(scan) == 1


def test_floor_scan_accepts_and_records_a_converged_solve():
    validator = _load_validator()
    solver = _Solver(converged=True, residual=7.0e-9, iters=311)
    scan = validator.floor_scan(solver, np.zeros(4), 1.0e3, [1.0], 1.0e-8, 1200)
    assert scan[0]["converged"] is True
    assert scan[0]["iterations"] == 311
    assert scan[0]["relative_residual"] == pytest.approx(7.0e-9)
    assert scan[0]["requested_tolerance"] == pytest.approx(1.0e-8)


def test_floor_scan_refuses_to_guess_when_the_solver_reports_no_convergence_data():
    validator = _load_validator()
    solver = _Solver(converged=True, residual=0.0)
    solver.timings = False
    with pytest.raises(validator.SolveContractError):
        validator.floor_scan(solver, np.zeros(4), 1.0e3, [1.0], 1.0e-8, 1200)


class _Spectrum:
    """A small generalized pair N v = lambda M v with M = I."""

    def __init__(self, eigenvalues):
        self.n = len(eigenvalues)
        self.N = np.diag(np.asarray(eigenvalues, float))

    def apply_configured_demag(self, x, _flag):
        return self.N @ np.asarray(x, float)

    def apply_configured_mass_riesz(self, x):
        return np.asarray(x, float)


def test_lobpcg_reports_non_convergence_instead_of_a_bare_ritz_value():
    validator = _load_validator()
    import scipy.sparse as sp

    n = 60
    eigenvalues = np.linspace(1.0e-9, 1.0, n)
    G = _Spectrum(eigenvalues)
    _, _, _, _, starved = validator.lobpcg_generalized(
        G, sp.identity(n, format="csr"), n, k=4, maxiter=1)
    assert starved["converged_low"] is False
    assert starved["final_residuals_low"]
    # SciPy's own non-convergence warning is captured, not swallowed.
    assert any("not reaching the requested tolerance" in w
               for w in starved["warnings_low"])

    _, _, _, _, settled = validator.lobpcg_generalized(
        G, sp.identity(n, format="csr"), n, k=4, maxiter=400)
    assert settled["tolerance"] == pytest.approx(1.0e-8)
    assert set(settled) >= {"converged_low", "converged_high",
                            "final_residuals_low", "final_residuals_high"}
