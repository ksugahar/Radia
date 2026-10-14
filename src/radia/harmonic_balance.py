"""Fail-closed utilities for hysteresis cycles and odd-harmonic balance.

The functions in this module are deliberately independent of a particular
field discretization.  A FEM, BEM, reduced-order, or analytical field update
can be supplied to :func:`solve_odd_harmonic_balance`; the sampling, Fourier
normalization, convergence history, and waveform-closure checks remain the
same.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Callable, Iterable

import numpy as np


Array = np.ndarray


@dataclass(frozen=True)
class HysteresisCycleMetrics:
    """Signed energy and closure diagnostics for one sampled B-H cycle."""

    signed_energy_density_j_per_m3: float
    closure_relative: float
    point_count: int
    passive_orientation: bool


@dataclass(frozen=True)
class OddHarmonicBalanceResult:
    """Result of a real, sine-phase, odd-harmonic fixed-point solve."""

    harmonics: tuple[int, ...]
    h_coefficients: Array
    b_coefficients: Array
    residual_history: Array
    iterations: int
    converged: bool
    waveform_closure_relative: float
    samples_per_period: int
    period_count: int


def _validated_harmonics(harmonics: Iterable[int]) -> tuple[int, ...]:
    values = tuple(int(value) for value in harmonics)
    if not values:
        raise ValueError("harmonics must not be empty")
    if values != tuple(sorted(set(values))):
        raise ValueError("harmonics must be unique and strictly increasing")
    if any(value <= 0 or value % 2 == 0 for value in values):
        raise ValueError("only positive odd harmonics are supported")
    return values


def hysteresis_cycle_metrics(
    field_a_per_m: Array,
    flux_density_t: Array,
    *,
    closure_scale_floor: float = 1.0e-30,
) -> HysteresisCycleMetrics:
    """Return signed ``integral H dB`` and endpoint closure for one cycle.

    The sign is not hidden with ``abs``.  A passive cycle sampled in the
    conventional time direction has positive ``integral H dB``.  Callers can
    therefore distinguish a valid loss from an accidentally reversed path.
    The supplied cycle should include both endpoints.
    """

    field = np.asarray(field_a_per_m, dtype=float)
    flux = np.asarray(flux_density_t, dtype=float)
    if field.ndim != 1 or flux.ndim != 1 or field.shape != flux.shape:
        raise ValueError("field and flux_density must be equal-length 1D arrays")
    if field.size < 4:
        raise ValueError("a cycle needs at least four samples")
    if not np.all(np.isfinite(field)) or not np.all(np.isfinite(flux)):
        raise ValueError("cycle samples must be finite")

    signed_energy = float(np.sum(0.5 * (field[:-1] + field[1:]) * np.diff(flux)))
    field_span = max(float(np.ptp(field)), closure_scale_floor)
    flux_span = max(float(np.ptp(flux)), closure_scale_floor)
    closure = float(
        np.hypot(
            (field[-1] - field[0]) / field_span,
            (flux[-1] - flux[0]) / flux_span,
        )
    )
    return HysteresisCycleMetrics(
        signed_energy_density_j_per_m3=signed_energy,
        closure_relative=closure,
        point_count=int(field.size),
        passive_orientation=bool(signed_energy > 0.0),
    )


def periodic_phase(samples_per_period: int, period_count: int = 1) -> Array:
    """Return endpoint-exclusive phase samples over an integer period count."""

    samples = int(samples_per_period)
    periods = int(period_count)
    if samples < 8:
        raise ValueError("samples_per_period must be at least 8")
    if periods < 1 or periods != period_count:
        raise ValueError("period_count must be a positive integer")
    count = samples * periods
    return 2.0 * np.pi * np.arange(count, dtype=float) / samples


def synthesize_odd_sine_series(
    coefficients: Array,
    harmonics: Iterable[int],
    phase: Array,
) -> Array:
    """Synthesize a real sine series with the requested odd harmonics."""

    modes = _validated_harmonics(harmonics)
    values = np.asarray(coefficients, dtype=float)
    angles = np.asarray(phase, dtype=float)
    if values.shape != (len(modes),):
        raise ValueError("coefficients must have one value per harmonic")
    if angles.ndim != 1 or not np.all(np.isfinite(angles)):
        raise ValueError("phase must be a finite 1D array")
    if not np.all(np.isfinite(values)):
        raise ValueError("coefficients must be finite")
    basis = np.sin(np.outer(np.asarray(modes, dtype=float), angles))
    return values @ basis


def project_odd_sine_harmonics(
    waveform: Array,
    harmonics: Iterable[int],
    *,
    samples_per_period: int,
    period_count: int = 1,
) -> Array:
    """Project endpoint-exclusive integer-period samples onto odd sine modes."""

    modes = _validated_harmonics(harmonics)
    phase = periodic_phase(samples_per_period, period_count)
    values = np.asarray(waveform, dtype=float)
    if values.shape != phase.shape:
        raise ValueError(
            "waveform length must equal samples_per_period * period_count; "
            "do not duplicate the final endpoint"
        )
    if not np.all(np.isfinite(values)):
        raise ValueError("waveform samples must be finite")
    basis = np.sin(np.outer(np.asarray(modes, dtype=float), phase))
    return (2.0 / values.size) * (basis @ values)


def solve_odd_harmonic_balance(
    material_response: Callable[[Array], Array],
    field_update: Callable[[Array], Array],
    initial_h_coefficients: Array,
    *,
    harmonics: Iterable[int] = (1, 3, 5),
    samples_per_period: int = 2048,
    period_count: int = 1,
    damping: float = 0.5,
    relative_tolerance: float = 1.0e-9,
    absolute_tolerance: float = 1.0e-12,
    max_iterations: int = 100,
    fail_on_nonconvergence: bool = True,
) -> OddHarmonicBalanceResult:
    """Solve a real odd-harmonic material/field fixed point.

    ``material_response`` maps the reconstructed H waveform to B samples.
    ``field_update`` maps the projected B coefficients to the next H
    coefficients.  This makes the routine usable with a full field solver
    without embedding solver-specific state here.  Only a converged iterate is
    returned; trial iterates remain local to this function.
    """

    modes = _validated_harmonics(harmonics)
    h_coefficients = np.asarray(initial_h_coefficients, dtype=float).copy()
    if h_coefficients.shape != (len(modes),) or not np.all(np.isfinite(h_coefficients)):
        raise ValueError("initial_h_coefficients must be finite and match harmonics")
    if not (0.0 < damping <= 1.0):
        raise ValueError("damping must lie in (0, 1]")
    if relative_tolerance < 0.0 or absolute_tolerance < 0.0:
        raise ValueError("tolerances must be non-negative")
    if int(max_iterations) < 1:
        raise ValueError("max_iterations must be positive")

    phase = periodic_phase(samples_per_period, period_count)
    residuals: list[float] = []
    converged = False

    for _ in range(int(max_iterations)):
        h_waveform = synthesize_odd_sine_series(h_coefficients, modes, phase)
        b_waveform = np.asarray(material_response(h_waveform), dtype=float)
        if b_waveform.shape != h_waveform.shape or not np.all(np.isfinite(b_waveform)):
            raise ValueError("material_response must return a finite waveform of matching shape")
        b_coefficients = project_odd_sine_harmonics(
            b_waveform,
            modes,
            samples_per_period=samples_per_period,
            period_count=period_count,
        )
        proposed = np.asarray(field_update(b_coefficients.copy()), dtype=float)
        if proposed.shape != h_coefficients.shape or not np.all(np.isfinite(proposed)):
            raise ValueError("field_update must return finite coefficients matching harmonics")
        next_coefficients = h_coefficients + damping * (proposed - h_coefficients)
        delta = float(np.linalg.norm(next_coefficients - h_coefficients))
        scale = max(float(np.linalg.norm(next_coefficients)), float(np.linalg.norm(h_coefficients)))
        residual = delta / max(scale, absolute_tolerance, np.finfo(float).tiny)
        residuals.append(residual)
        h_coefficients = next_coefficients
        if delta <= absolute_tolerance + relative_tolerance * scale:
            converged = True
            break

    h_waveform = synthesize_odd_sine_series(h_coefficients, modes, phase)
    b_waveform = np.asarray(material_response(h_waveform), dtype=float)
    b_coefficients = project_odd_sine_harmonics(
        b_waveform,
        modes,
        samples_per_period=samples_per_period,
        period_count=period_count,
    )
    reconstructed_b = synthesize_odd_sine_series(b_coefficients, modes, phase)
    waveform_scale = max(float(np.linalg.norm(b_waveform)), np.finfo(float).tiny)
    waveform_closure = float(np.linalg.norm(b_waveform - reconstructed_b) / waveform_scale)

    result = OddHarmonicBalanceResult(
        harmonics=modes,
        h_coefficients=h_coefficients.copy(),
        b_coefficients=b_coefficients.copy(),
        residual_history=np.asarray(residuals, dtype=float),
        iterations=len(residuals),
        converged=converged,
        waveform_closure_relative=waveform_closure,
        samples_per_period=int(samples_per_period),
        period_count=int(period_count),
    )
    if fail_on_nonconvergence and not converged:
        final = residuals[-1] if residuals else float("inf")
        raise RuntimeError(
            f"odd-harmonic balance did not converge in {max_iterations} iterations "
            f"(relative coefficient residual {final:.3e})"
        )
    return result


__all__ = [
    "HysteresisCycleMetrics",
    "OddHarmonicBalanceResult",
    "hysteresis_cycle_metrics",
    "periodic_phase",
    "project_odd_sine_harmonics",
    "solve_odd_harmonic_balance",
    "synthesize_odd_sine_series",
]
