"""Fail-loud contracts for FEM-Kelvin solver mode and Karl acceptance."""

import sys
from pathlib import Path

import pytest

PANELS = Path(__file__).parents[1] / "src" / "radia" / "panels"
sys.path.insert(0, str(PANELS))
import calc_fem_kelvin as solver  # noqa: E402


def _validate(**overrides):
    arguments = dict(
        solver="sparsecholesky", formulation="total", max_iter=4,
        esim_per_panel=False, periodic_kelvin=False)
    arguments.update(overrides)
    solver._validate_mode_contracts(**arguments)


@pytest.mark.parametrize("overrides, message", [
    ({"max_iter": 0}, "positive"),
    ({"esim_per_panel": True}, "local-Z surface loss"),
    ({"solver": "iccg", "periodic_kelvin": True}, "periodic Kelvin"),
    ({"solver": "bddc", "formulation": "scattered"}, "only for.*sparsecholesky"),
])
def test_unsupported_modes_fail_before_solving(overrides, message):
    with pytest.raises(ValueError, match=message):
        _validate(**overrides)


def test_supported_total_and_direct_scattered_modes_remain_available():
    _validate(solver="bddc", formulation="total")
    _validate(solver="sparsecholesky", formulation="scattered",
              periodic_kelvin=True)


def test_karl_acceptance_publishes_impedance_used_by_final_field():
    converged, published = solver._accepted_karl_impedance(
        iteration=2, max_iter=5, relative_change=1.0e-4,
        tolerance=1.0e-3, impedance_used_for_solve=2.0 + 3.0j)
    assert converged
    assert published == 2.0 + 3.0j

    converged, published = solver._accepted_karl_impedance(
        iteration=4, max_iter=5, relative_change=2.0e-3,
        tolerance=1.0e-3, impedance_used_for_solve=2.0 + 3.0j)
    assert not converged
    assert published is None


def test_fem_cell_convergence_gate_keeps_diagnostics():
    with pytest.raises(
            RuntimeError,
            match=r"FEM scalar ESIM.*iterations=9.*relative_change=0.02"):
        solver.require_esim_converged(
            {"converged": False, "iterations": 9,
             "relative_change": 0.02, "Z": 1.0 + 1.0j},
            "FEM scalar ESIM cell solve")
