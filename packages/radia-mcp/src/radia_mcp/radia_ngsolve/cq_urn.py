"""Convolution-quadrature bridge for URN relaxation models.

The Universal Relaxation Network (URN) identifies passive, causal relaxation
terms from frequency-domain data.  Convolution quadrature (CQ) then needs only
the Laplace-domain evaluator ``H(s)`` of that fitted model.  These helpers keep
that bridge explicit and small enough to use as teaching material before it is
embedded in a full FEM/BEM or Maxwell time-domain solver.
"""

from __future__ import annotations

from dataclasses import dataclass
import datetime as _datetime
import importlib.metadata
import math
import platform
import sys
import time
from collections.abc import Callable, Sequence

import numpy as np


def _real_1d(name: str, values: Sequence[float]) -> np.ndarray:
    arr = np.asarray(values, dtype=float).ravel()
    if arr.size == 0:
        raise ValueError(f"{name} must contain at least one value")
    if not np.all(np.isfinite(arr)):
        raise ValueError(f"{name} must be finite")
    return arr


def relaxation_response(
    s: np.ndarray | complex | float,
    weights: Sequence[float],
    taus: Sequence[float],
    feedthrough: float = 0.0,
) -> np.ndarray:
    """Evaluate a passive Debye-relaxation surrogate.

    ``H(s) = feedthrough + sum_k weights[k] / (1 + s*tau[k])``.

    This is the simplest URN-exportable model.  General URN components such as
    Cole-Cole or CPE terms can be approximated by a short positive Debye ladder,
    which is the same structure used by ADE/SPICE realizations.
    """

    w = _real_1d("weights", weights)
    tau = _real_1d("taus", taus)
    if w.shape != tau.shape:
        raise ValueError("weights and taus must have the same length")
    if np.any(w < 0.0):
        raise ValueError("weights must be non-negative for the passive bridge")
    if np.any(tau <= 0.0):
        raise ValueError("taus must be positive")
    z = np.asarray(s, dtype=complex)
    out = np.zeros_like(z, dtype=complex) + float(feedthrough)
    for wk, tk in zip(w, tau, strict=True):
        out = out + wk / (1.0 + z * tk)
    return out


@dataclass(frozen=True)
class NonnegativeDebyeFit:
    """A compact URN-compatible non-negative Debye ladder fit."""

    feedthrough: float
    weights: tuple[float, ...]
    taus: tuple[float, ...]
    relative_error: float
    active_count: int
    active_threshold: float

    def evaluate(self, s: np.ndarray | complex | float) -> np.ndarray:
        return relaxation_response(s, self.weights, self.taus, self.feedthrough)

    def as_dict(self) -> dict:
        return {
            "feedthrough": self.feedthrough,
            "weights": list(self.weights),
            "taus_s": list(self.taus),
            "relative_error": self.relative_error,
            "active_count": self.active_count,
            "active_threshold": self.active_threshold,
        }


def fit_nonnegative_debye(
    freq_hz: Sequence[float],
    response: Sequence[complex],
    tau_grid: Sequence[float],
    *,
    fit_feedthrough: bool = True,
    active_threshold: float = 1.0e-8,
) -> NonnegativeDebyeFit:
    """Fit ``response(j*omega)`` by a non-negative Debye ladder.

    The routine is intentionally modest: it is not a replacement for the full
    URN trainer, but it exposes the post-fit contract that CQ consumes.  When
    SciPy is available it uses NNLS; otherwise it falls back to a clipped least
    squares solution so the teaching example remains runnable.
    """

    freq = _real_1d("freq_hz", freq_hz)
    tau = _real_1d("tau_grid", tau_grid)
    if np.any(freq <= 0.0):
        raise ValueError("freq_hz must be positive")
    if np.any(tau <= 0.0):
        raise ValueError("tau_grid must be positive")
    target = np.asarray(response, dtype=complex).ravel()
    if target.shape != freq.shape:
        raise ValueError("freq_hz and response must have the same length")

    s = 1j * 2.0 * np.pi * freq
    columns: list[np.ndarray] = []
    if fit_feedthrough:
        columns.append(np.ones_like(target, dtype=complex))
    for tk in tau:
        columns.append(1.0 / (1.0 + s * tk))
    basis = np.column_stack(columns)
    a = np.vstack([basis.real, basis.imag])
    b = np.concatenate([target.real, target.imag])

    coeff: np.ndarray
    try:
        from scipy.optimize import nnls  # type: ignore

        coeff, _ = nnls(a, b)
    except Exception:  # pragma: no cover - exercised only without SciPy/NNLS
        coeff = np.linalg.lstsq(a, b, rcond=None)[0]
        coeff = np.maximum(coeff, 0.0)

    if fit_feedthrough:
        feedthrough = float(coeff[0])
        weights = coeff[1:]
    else:
        feedthrough = 0.0
        weights = coeff

    fitted = relaxation_response(s, weights, tau, feedthrough)
    denom = np.linalg.norm(target)
    rel = float(np.linalg.norm(fitted - target) / (denom + 1.0e-30))
    active = int(np.count_nonzero(np.asarray(weights) > active_threshold))
    if feedthrough > active_threshold:
        active += 1
    return NonnegativeDebyeFit(
        feedthrough=feedthrough,
        weights=tuple(float(x) for x in weights),
        taus=tuple(float(x) for x in tau),
        relative_error=rel,
        active_count=active,
        active_threshold=float(active_threshold),
    )


def bdf_delta(zeta: np.ndarray | complex | float, method: str = "bdf2") -> np.ndarray:
    """Return the BDF generating polynomial ``delta(zeta)`` used by CQ."""

    z = np.asarray(zeta, dtype=complex)
    m = method.strip().lower()
    if m in {"bdf1", "backward_euler", "be"}:
        return 1.0 - z
    if m == "bdf2":
        return 1.5 - 2.0 * z + 0.5 * z * z
    raise ValueError("method must be 'bdf1' or 'bdf2'")


def next_power_of_two(n: int) -> int:
    """Return the smallest power of two greater than or equal to ``n``."""

    if n <= 1:
        return 1
    return 1 << (int(n) - 1).bit_length()


def cq_time_grid_contract_gate(
    dt: float,
    n_steps: int,
    *,
    method: str = "bdf2",
    contour_samples: int | None = None,
    radius: float | None = None,
) -> dict:
    """Validate the implementation-neutral Lubich CQ time-grid contract."""

    dt_value = float(dt)
    n_value = int(n_steps)
    m_value = int(contour_samples) if contour_samples is not None else next_power_of_two(2 * n_value)
    rho = float(radius) if radius is not None else np.finfo(float).eps ** (0.5 / m_value)
    method_key = str(method or "").strip().lower()
    method_ok = method_key in {"bdf1", "bdf2"}
    shape_ok = dt_value > 0.0 and n_value >= 2 and m_value >= n_value and 0.0 < rho < 1.0
    min_real_s = float("nan")
    if method_ok and shape_ok:
        j = np.arange(m_value, dtype=float)
        zeta = rho * np.exp(-2j * np.pi * j / m_value)
        min_real_s = float(np.min(np.real(bdf_delta(zeta, method=method_key) / dt_value)))
    checks = {
        "method_supported": method_ok,
        "time_step_positive": dt_value > 0.0,
        "step_count_valid": n_value >= 2,
        "contour_covers_time_grid": m_value >= n_value,
        "radius_inside_unit_circle": 0.0 < rho < 1.0,
        "laplace_nodes_in_right_half_plane": np.isfinite(min_real_s) and min_real_s > 0.0,
    }
    return {
        "policy": "cq_time_grid_contract_gate",
        "status": "ok" if all(checks.values()) else "needs_attention",
        "method": method_key,
        "time_step": dt_value,
        "num_time_steps": n_value,
        "contour_samples": m_value,
        "radius": rho,
        "time_end": (n_value - 1) * dt_value,
        "min_real_laplace_node": min_real_s,
        "checks": checks,
        "notes": [
            "M>=N is shared; M=N and padded power-of-two FFTs are both valid implementations.",
            "Compare MATLAB and Python CQ on physical observables, not identical FFT padding choices.",
        ],
    }


def cq_weights_from_laplace(
    laplace_response: Callable[[np.ndarray], np.ndarray],
    dt: float,
    n_steps: int,
    *,
    method: str = "bdf2",
    radius: float | None = None,
    fft_size: int | None = None,
) -> np.ndarray:
    """Compute Lubich CQ weights from a Laplace-domain response ``H(s)``.

    The caller supplies an evaluator for the already-identified passive model.
    No time-domain state realization is assumed here, which makes the same code
    usable for acoustic FEM/BEM, impedance boundaries, and Maxwell material
    kernels.
    """

    if dt <= 0.0:
        raise ValueError("dt must be positive")
    if n_steps <= 0:
        raise ValueError("n_steps must be positive")
    n_fft = int(fft_size) if fft_size is not None else next_power_of_two(2 * n_steps)
    if n_fft < n_steps:
        raise ValueError("fft_size must be >= n_steps")
    rho = float(radius) if radius is not None else np.finfo(float).eps ** (0.5 / n_fft)
    if not (0.0 < rho < 1.0):
        raise ValueError("radius must be in (0, 1)")

    j = np.arange(n_fft, dtype=float)
    zeta = rho * np.exp(-2j * np.pi * j / n_fft)
    s = bdf_delta(zeta, method=method) / float(dt)
    values = np.asarray(laplace_response(s), dtype=complex)
    if values.shape != s.shape:
        values = np.broadcast_to(values, s.shape)
    coeff = np.fft.ifft(values)
    scale = rho ** np.arange(n_steps)
    return coeff[:n_steps] / scale


def cq_convolve(weights: Sequence[complex], input_signal: Sequence[float]) -> np.ndarray:
    """Apply CQ weights to a sampled input by causal discrete convolution."""

    w = np.asarray(weights, dtype=complex).ravel()
    u = np.asarray(input_signal, dtype=float).ravel()
    if w.size == 0 or u.size == 0:
        raise ValueError("weights and input_signal must be non-empty")
    return np.convolve(w, u)[: u.size]


def periodic_ifft_response(
    laplace_response: Callable[[np.ndarray], np.ndarray],
    dt: float,
    input_signal: Sequence[float],
) -> np.ndarray:
    """Naive periodic frequency-domain response used only as a teaching contrast."""

    u = np.asarray(input_signal, dtype=float).ravel()
    omega = 2.0 * np.pi * np.fft.fftfreq(u.size, d=dt)
    values = np.asarray(laplace_response(1j * omega), dtype=complex)
    return np.fft.ifft(values * np.fft.fft(u))


def _version_info() -> dict:
    versions = {
        "python": sys.version.split()[0],
        "platform": platform.platform(),
    }
    for dist_name, key in [
        ("radia", "radia_version"),
        ("radia-mcp", "radia_mcp_version"),
        ("numpy", "numpy_version"),
        ("scipy", "scipy_version"),
    ]:
        try:
            versions[key] = importlib.metadata.version(dist_name)
        except Exception:
            versions[key] = None
    return versions


def make_cq_urn_bridge_artifact(
    *,
    n_steps: int = 100,
    dt: float = 0.01,
    hit_index: int = 10,
    cq_method: str = "bdf2",
) -> dict:
    """Build a compact URN -> CQ teaching artifact with numerical checks."""

    start = time.perf_counter()
    if hit_index < 0 or hit_index >= n_steps:
        raise ValueError("hit_index must lie within the sampled time window")

    true_weights = np.array([0.62, 0.30])
    true_taus = np.array([0.02, 0.35])
    feedthrough = 0.08

    freq_hz = np.logspace(-1.0, 2.0, 80)
    target = relaxation_response(
        1j * 2.0 * np.pi * freq_hz,
        true_weights,
        true_taus,
        feedthrough,
    )
    tau_grid = np.array([0.005, 0.01, 0.02, 0.05, 0.1, 0.35, 0.7, 1.5])
    t_fit0 = time.perf_counter()
    fit = fit_nonnegative_debye(freq_hz, target, tau_grid, active_threshold=1.0e-6)
    fit_s = time.perf_counter() - t_fit0

    laplace = fit.evaluate
    t_cq0 = time.perf_counter()
    weights = cq_weights_from_laplace(laplace, dt, n_steps, method=cq_method)
    cq_weights_s = time.perf_counter() - t_cq0

    time_s = np.arange(n_steps, dtype=float) * dt
    input_signal = np.zeros(n_steps, dtype=float)
    input_signal[hit_index:] = 1.0
    t_conv0 = time.perf_counter()
    cq_response = cq_convolve(weights, input_signal).real
    ifft_response = periodic_ifft_response(laplace, dt, input_signal).real
    convolution_s = time.perf_counter() - t_conv0

    omega = 2.0 * np.pi * np.logspace(-1.0, 3.0, 128)
    freq_axis_real_min = float(np.min(laplace(1j * omega).real))
    positive_axis_real_min = float(np.min(laplace(np.logspace(-1.0, 3.0, 128)).real))
    prehit_cq = float(np.max(np.abs(cq_response[:hit_index]))) if hit_index else 0.0
    prehit_ifft = float(np.max(np.abs(ifft_response[:hit_index]))) if hit_index else 0.0
    after = slice(hit_index, n_steps)
    rms_scale = math.sqrt(float(np.mean(cq_response[after] ** 2))) + 1.0e-30
    rms_ifft_vs_cq = float(
        math.sqrt(float(np.mean((ifft_response[after] - cq_response[after]) ** 2)))
        / rms_scale
    )

    total_s = time.perf_counter() - start
    artifact = {
        "schema": "radia.cq_urn_bridge.v1",
        "result_artifact_id": "docs.universal_relaxation_network.cq_urn_bridge",
        "generated_at_utc": _datetime.datetime.now(
            _datetime.timezone.utc
        ).isoformat(timespec="seconds").replace("+00:00", "Z"),
        "versions": _version_info(),
        "lesson": (
            "A passive URN relaxation fit can be used directly as H(s) in "
            "Lubich convolution quadrature.  CQ gives a causal time-domain "
            "operator without imposing the periodic wrap-around implicit in a "
            "plain FFT/IFFT demonstration."
        ),
        "model": {
            "type": "nonnegative_debye_ladder",
            "true_feedthrough": feedthrough,
            "true_weights": true_weights.tolist(),
            "true_taus_s": true_taus.tolist(),
            "candidate_taus_s": tau_grid.tolist(),
            "fit": fit.as_dict(),
        },
        "cq": {
            "method": cq_method,
            "dt_s": dt,
            "n_steps": n_steps,
            "hit_index": hit_index,
            "hit_time_s": hit_index * dt,
            "prehit_max_abs": prehit_cq,
        },
        "ifft_periodic_contrast": {
            "prehit_max_abs": prehit_ifft,
            "rms_difference_after_hit_relative_to_cq": rms_ifft_vs_cq,
        },
        "checks": {
            "fit_relative_error_lt_1e-10": bool(fit.relative_error < 1.0e-10),
            "nonnegative_fit_weights": bool(np.all(np.asarray(fit.weights) >= -1.0e-14)),
            "positive_real_frequency_axis_min": freq_axis_real_min,
            "positive_real_positive_axis_min": positive_axis_real_min,
            "cq_causal_before_hit": bool(prehit_cq < 1.0e-12),
            "ifft_periodic_has_prehit_wraparound": bool(prehit_ifft > 1.0e-3),
        },
        "result_output_schema_id": "radia.cq_urn_bridge.timeseries.v1",
        "result_output_columns": [
            "time_s",
            "input_unit_step",
            "cq_response",
            "ifft_periodic_response",
        ],
        "result_output_units": {
            "time_s": "s",
            "input_unit_step": "1",
            "cq_response": "normalized",
            "ifft_periodic_response": "normalized",
        },
        "timeseries": {
            "time_s": time_s.tolist(),
            "input_unit_step": input_signal.tolist(),
            "cq_response": cq_response.tolist(),
            "ifft_periodic_response": ifft_response.tolist(),
        },
        "timing_breakdown_s": {
            "fit_nonnegative_debye": fit_s,
            "cq_weight_generation": cq_weights_s,
            "convolution_and_ifft_contrast": convolution_s,
            "total": total_s,
        },
        "mcp_feedback": {
            "public_summary": (
                "CQ consumes the URN fit through a Laplace evaluator H(s); "
                "for time-domain acoustic or Maxwell demos, prefer CQ for "
                "causal transients and keep IFFT as a periodic steady-state "
                "contrast unless zero-padding/windowing is explicitly handled."
            ),
            "learning_targets": ["radia-mcp:radia_ngsolve.urn(topic='cq')"],
            "learning_lanes": {"public": "encoded", "source_tool": "none"},
        },
    }
    artifact["pass"] = bool(
        artifact["checks"]["fit_relative_error_lt_1e-10"]
        and artifact["checks"]["nonnegative_fit_weights"]
        and artifact["checks"]["cq_causal_before_hit"]
        and artifact["checks"]["ifft_periodic_has_prehit_wraparound"]
    )
    return artifact
