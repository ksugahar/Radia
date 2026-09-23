"""Contracts for the terminal-driven A-V eddy solver.

The DC limit is the sharpest cheap test this formulation has: as the frequency
falls the skin effect vanishes, A stops influencing the current, and the
terminal resistance must collapse onto ``L / (sigma A)`` of the meshed
conductor with nothing left to tune.  Every bug found while building the
solver showed up here first -- a terminal face picked by centre-of-face radius,
a cap-flux integral that evaluated ``grad V`` to exactly zero on the boundary,
a symmetrising substitution that was never undone, and a vector potential left
with the natural curl-curl condition on the terminal caps, which cost 19%.

These run on a deliberately coarse mesh: the point is the contract, not the
accuracy, and the accuracy is measured against the Bessel solution in
``validation_test/induction_heating/aphi_round_wire_validation.py``.
"""
from __future__ import annotations

import math
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "validation_test" / "induction_heating"))

MU0 = 4e-7 * math.pi
RADIUS_M = 1.0e-3
LENGTH_M = 6.0e-3
SIGMA = 5.8e7


@pytest.fixture(scope="module")
def wire():
    pytest.importorskip("ngsolve")
    pytest.importorskip("netgen.occ")
    from aphi_round_wire_validation import build_wire_mesh

    return build_wire_mesh(RADIUS_M, LENGTH_M, 4.0, 0.6e-3, 3.0e-3, 0, 0.5)


def _dc_resistance(mesh):
    import ngsolve as ng

    volume = ng.Integrate(ng.CF(1.0), mesh,
                          definedon=mesh.Materials("conductor")).real
    return LENGTH_M * LENGTH_M / (SIGMA * volume)


def _solve(mesh, a_over_delta, **kwargs):
    from radia.eddy_aphi import solve_eddy_aphi

    delta = RADIUS_M / a_over_delta
    omega = 2.0 / (delta * delta * MU0 * SIGMA)
    return solve_eddy_aphi(mesh, conductor="conductor", source="source",
                           sink="sink", frequency_hz=omega / (2.0 * math.pi),
                           sigma=SIGMA, **kwargs)


def test_dc_limit_is_the_meshed_conductor_resistance(wire):
    result = _solve(wire, 1e-3, order=2)
    assert result.resistance_ohm == pytest.approx(_dc_resistance(wire),
                                                  rel=2e-3)


@pytest.mark.parametrize("order", [1, 2])
def test_dc_limit_holds_at_every_order(wire, order):
    result = _solve(wire, 1e-3, order=order)
    assert result.resistance_ohm == pytest.approx(_dc_resistance(wire),
                                                  rel=5e-3)


def test_resistance_rises_with_frequency(wire):
    """Skin effect, in the one direction it can possibly go."""
    low = _solve(wire, 1e-2, order=2).resistance_ohm
    high = _solve(wire, 3.0, order=2).resistance_ohm
    assert high > low * 1.3


def test_terminal_current_is_normalised_to_one_amp(wire):
    """Every reported quantity is rescaled, so the scale must undo the drive."""
    result = _solve(wire, 1.0, order=2)
    assert abs(result.terminal_current_A * result.scale) == pytest.approx(
        1.0, rel=1e-12)
    assert result.terminal_impedance_ohm == pytest.approx(
        complex(result.resistance_ohm,
                result.inductance_H * 2.0 * math.pi * result.frequency_hz),
        rel=1e-12)


def test_zero_and_negative_frequency_are_refused(wire):
    from radia.eddy_aphi import solve_eddy_aphi

    for bad in (0.0, -1.0):
        with pytest.raises(ValueError, match="frequency_hz must be positive"):
            solve_eddy_aphi(wire, conductor="conductor", source="source",
                            sink="sink", frequency_hz=bad, sigma=SIGMA)


def test_unknown_region_and_boundary_names_are_refused(wire):
    from radia.eddy_aphi import solve_eddy_aphi

    with pytest.raises(ValueError, match="material name"):
        solve_eddy_aphi(wire, conductor="not_a_material", source="source",
                        sink="sink", frequency_hz=1e3, sigma=SIGMA)
    with pytest.raises(ValueError, match="boundary carries the name"):
        solve_eddy_aphi(wire, conductor="conductor", source="not_a_boundary",
                        sink="sink", frequency_hz=1e3, sigma=SIGMA)
