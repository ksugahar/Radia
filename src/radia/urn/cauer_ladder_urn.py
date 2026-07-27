"""Differentiable Cauer-ladder URN candidate.

The SA/RM-2026 Y-domain URN sums basis functions in parallel.  This module
explores the complementary direction discussed during review: keep the model
small, but increase expressiveness through a Cauer-like series/shunt continued
fraction.  The topology is fixed and all element values are positive, so the
network is directly differentiable and passive by construction.
"""

from __future__ import annotations

from dataclasses import dataclass, replace
from typing import Any

import numpy as np
import torch
from torch import nn, optim

from .y_admittance_urn import (
    _EPS,
    _as_1d_complex128,
    _as_1d_float64,
    _inv_sigmoid,
    _positive_range_from_raw,
    complex_smooth_l1,
    scattering_transform,
)


@dataclass
class CauerLadderURNConfig:
    """Configuration for a positive differentiable Cauer ladder.

    The model uses normalized impedance/admittance variables:

    ``zbar = Z / z_ref`` and ``x = s / omega_ref``.

    Each section is

    ``Z_k = r_k + l_k x + 1 / (g_k + c_k x + 1 / Z_{k+1})``.

    Thus a 5-section ladder has 20 positive parameters and a 6-section ladder
    has 24, close to the original 22-basis research-meeting dictionary.
    """

    n_sections: int = 6
    lr: float = 3.0e-3
    n_epochs: int = 3000
    n_restarts: int = 1
    seed: int = 17

    huber_delta: float = 1.0e-2
    regularization_weight: float = 1.0e-7
    log_smoothness_weight: float = 1.0e-7

    parameter_min: float = 1.0e-8
    parameter_max: float = 1.0e8
    init_series_resistance: float = 0.05
    init_series_inductance: float = 0.05
    init_shunt_conductance: float = 0.05
    init_shunt_capacitance: float = 0.05
    use_peeling_initialization: bool = False
    use_rational_initialization: bool = False
    rational_order: int = 6
    rational_sk_iterations: int = 4
    rational_regularization: float = 1.0e-10
    rational_cauer_max_nfev: int = 1200
    use_least_squares_polish: bool = False
    least_squares_max_nfev: int = 1200
    least_squares_z_weight: float = 1.0
    least_squares_y_weight: float = 1.0
    extra_series_impedance: float = 1.0e6
    extra_shunt_admittance: float = 1.0e-6

    z0: float | None = None
    z_ref: float | None = None
    omega_ref: float | None = None

    progressive_start_sections: int = 1
    progressive_stage_epochs: int | None = None
    progressive_lr_decay: float = 0.85

    alternating_cycles: int = 30
    alternating_block_epochs: int = 25
    alternating_reject_growth: float = 1.25
    alternating_lr_decay_on_reject: float = 0.5
    alternating_min_lr: float = 1.0e-6
    frozen_outer_sections: int = 0
    tail_train_cycles: int | None = None
    polish_train_cycles: int | None = None

    @property
    def total_parameters(self) -> int:
        return 4 * self.n_sections

    @classmethod
    def twenty_two_parameter_candidate(cls, **overrides: Any) -> "CauerLadderURNConfig":
        """Return a compact Cauer candidate with 24 positive parameters."""

        values = {
            "n_sections": 6,
            "n_epochs": 2500,
            "lr": 3.0e-3,
            "n_restarts": 1,
            "regularization_weight": 1.0e-7,
            "log_smoothness_weight": 1.0e-7,
        }
        values.update(overrides)
        return cls(**values)


def _raw_for_positive_value(value: float, lo: float, hi: float) -> float:
    value = float(np.clip(value, lo * (1.0 + 1.0e-9), hi * (1.0 - 1.0e-9)))
    frac = (np.log(value) - np.log(lo)) / (np.log(hi) - np.log(lo))
    return float(_inv_sigmoid(frac))


def _initial_raw(value: float, n: int, lo: float, hi: float) -> torch.Tensor:
    return torch.full((n,), _raw_for_positive_value(value, lo, hi), dtype=torch.float64)


def _polyval_ascending(coeffs: np.ndarray, x: np.ndarray) -> np.ndarray:
    out = np.zeros_like(x, dtype=np.complex128)
    for coeff in reversed(np.asarray(coeffs, dtype=np.complex128)):
        out = out * x + coeff
    return out


@dataclass(frozen=True)
class CauerRationalFit:
    """Small pole/zero rational fit used as a Cauer initializer."""

    numerator: np.ndarray
    denominator: np.ndarray
    z_ref: float
    omega_ref: float

    def predict(self, freqs_hz: np.ndarray | list[float]) -> np.ndarray:
        freqs = _as_1d_float64(freqs_hz, "freqs_hz")
        x = 1j * (2.0 * np.pi * freqs) / self.omega_ref
        den = _polyval_ascending(self.denominator, x)
        num = _polyval_ascending(self.numerator, x)
        return self.z_ref * num / (den + _EPS)

    @property
    def poles(self) -> np.ndarray:
        return np.roots(np.asarray(self.denominator, dtype=np.complex128)[::-1])

    @property
    def zeros(self) -> np.ndarray:
        return np.roots(np.asarray(self.numerator, dtype=np.complex128)[::-1])


def fit_rational_pole_zero(
    freqs_hz: np.ndarray | list[float],
    z_data: np.ndarray | list[complex],
    *,
    order: int = 6,
    omega_ref: float | None = None,
    z_ref: float | None = None,
    sk_iterations: int = 4,
    regularization: float = 1.0e-10,
) -> CauerRationalFit:
    """Fit ``Z/z_ref = P(x)/Q(x)`` with a small SK-style least-squares loop."""

    freqs = _as_1d_float64(freqs_hz, "freqs_hz")
    z_arr = _as_1d_complex128(z_data, "z_data")
    if freqs.shape != z_arr.shape:
        raise ValueError("freqs_hz and z_data must have the same length")
    if order < 1:
        raise ValueError("order must be positive")

    omega = 2.0 * np.pi * freqs
    w_ref = float(omega_ref or np.sqrt(float(np.min(omega)) * float(np.max(omega))))
    z_scale = float(z_ref or np.median(np.abs(z_arr)))
    w_ref = max(w_ref, _EPS)
    z_scale = max(z_scale, _EPS)

    x = 1j * omega / w_ref
    zbar = z_arr / z_scale
    num_degree = order
    den_degree = order
    num_vander = np.column_stack([x**k for k in range(num_degree + 1)])
    den_vander = np.column_stack([x**k for k in range(1, den_degree + 1)])
    den = np.ones_like(zbar)
    sol = np.zeros(num_degree + den_degree + 1, dtype=np.complex128)

    base_weight = 1.0 / np.maximum(np.abs(zbar), np.median(np.abs(zbar)) * 1.0e-6)
    for _ in range(max(1, int(sk_iterations))):
        a = np.hstack([num_vander, -zbar[:, None] * den_vander])
        b = zbar
        weight = base_weight / np.maximum(np.abs(den), 1.0e-12)
        aw = a * weight[:, None]
        bw = b * weight
        if regularization > 0.0:
            reg = np.sqrt(float(regularization)) * np.eye(a.shape[1], dtype=np.complex128)
            aw = np.vstack([aw, reg])
            bw = np.concatenate([bw, np.zeros(a.shape[1], dtype=np.complex128)])
        sol, *_ = np.linalg.lstsq(aw, bw, rcond=None)
        denominator = np.concatenate(
            [np.ones(1, dtype=np.complex128), sol[num_degree + 1 :]]
        )
        den = _polyval_ascending(denominator, x)

    numerator = sol[: num_degree + 1]
    denominator = np.concatenate([np.ones(1, dtype=np.complex128), sol[num_degree + 1 :]])
    return CauerRationalFit(
        numerator=np.asarray(numerator, dtype=np.complex128),
        denominator=np.asarray(denominator, dtype=np.complex128),
        z_ref=z_scale,
        omega_ref=w_ref,
    )


def _fit_positive_rx(values: np.ndarray, x_imag: np.ndarray) -> tuple[float, float]:
    finite = np.isfinite(values.real) & np.isfinite(values.imag)
    if not np.any(finite):
        return 0.0, 0.0
    real = values.real[finite]
    imag = values.imag[finite]
    x = x_imag[finite]
    positive_real = real[real > 0.0]
    r = float(np.percentile(positive_real, 10)) if positive_real.size else 0.0
    slope_candidates = imag[x > 0.0] / x[x > 0.0]
    positive_slope = slope_candidates[np.isfinite(slope_candidates) & (slope_candidates > 0.0)]
    if positive_slope.size:
        l = float(np.percentile(positive_slope, 10))
    else:
        denom = float(np.dot(x, x))
        l = float(np.dot(x, imag) / denom) if denom > 0.0 else 0.0
    return max(r, 0.0), max(l, 0.0)


def _stable_subtract(
    values: np.ndarray,
    element: np.ndarray,
    *,
    min_relative_residual: float = 1.0e-4,
) -> tuple[float, np.ndarray]:
    scale_candidates = (1.0, 0.5, 0.25, 0.125, 0.0625, 0.0)
    reference = max(float(np.median(np.abs(values))), 1.0e-24)
    for scale in scale_candidates:
        residual = values - scale * element
        if not np.all(np.isfinite(residual.real)) or not np.all(np.isfinite(residual.imag)):
            continue
        if float(np.percentile(np.abs(residual), 10)) < min_relative_residual * reference:
            continue
        return float(scale), residual
    return 0.0, values


def _cauer_numpy_from_positive_params(
    params: np.ndarray,
    x: np.ndarray,
    n_sections: int,
) -> np.ndarray:
    r = params[0:n_sections]
    l = params[n_sections : 2 * n_sections]
    g = params[2 * n_sections : 3 * n_sections]
    c = params[3 * n_sections : 4 * n_sections]
    z_next: np.ndarray | None = None
    for i in range(n_sections - 1, -1, -1):
        y_tail = 0.0 if z_next is None else 1.0 / (z_next + _EPS)
        y_shunt = g[i] + c[i] * x + y_tail
        z_next = r[i] + l[i] * x + 1.0 / (y_shunt + _EPS)
    if z_next is None:
        raise RuntimeError("empty Cauer ladder")
    return np.asarray(z_next, dtype=np.complex128)


def _relative_complex_residual(pred: np.ndarray, target: np.ndarray) -> np.ndarray:
    floor = max(float(np.median(np.abs(target))) * 1.0e-6, 1.0e-24)
    scale = np.maximum(np.abs(target), floor)
    residual = (pred - target) / scale
    return np.concatenate([residual.real, residual.imag])


class CauerLadderURN(nn.Module):
    """Positive Cauer ladder evaluated by differentiable continued fractions."""

    def __init__(
        self,
        freqs_hz: np.ndarray | list[float],
        config: CauerLadderURNConfig | None = None,
        *,
        z_data: np.ndarray | list[complex] | None = None,
    ) -> None:
        super().__init__()
        self.config = config or CauerLadderURNConfig()
        if self.config.n_sections <= 0:
            raise ValueError("n_sections must be positive")
        if self.config.parameter_min <= 0.0 or self.config.parameter_max <= self.config.parameter_min:
            raise ValueError("parameter_min/parameter_max must be positive and ordered")

        freqs = _as_1d_float64(freqs_hz, "freqs_hz")
        if np.any(freqs <= 0.0):
            raise ValueError("all frequencies must be positive")
        self.freqs_hz = freqs
        self.omega = 2.0 * np.pi * freqs
        self.omega_min = float(np.min(self.omega))
        self.omega_max = float(np.max(self.omega))
        self.omega_ref = float(
            self.config.omega_ref or np.sqrt(self.omega_min * self.omega_max)
        )

        if z_data is not None:
            z_arr = _as_1d_complex128(z_data, "z_data")
            if z_arr.shape != freqs.shape:
                raise ValueError("freqs_hz and z_data must have the same length")
            z_ref = self.config.z_ref or float(np.median(np.abs(z_arr)))
            z0 = self.config.z0 or float(np.median(np.abs(z_arr)))
        else:
            z_ref = self.config.z_ref or 1.0
            z0 = self.config.z0 or 1.0
        self.z_ref = float(max(z_ref, _EPS))
        self.z0 = float(max(z0, _EPS))

        n = self.config.n_sections
        lo = self.config.parameter_min
        hi = self.config.parameter_max
        self.r_raw = nn.Parameter(_initial_raw(self.config.init_series_resistance, n, lo, hi))
        self.l_raw = nn.Parameter(_initial_raw(self.config.init_series_inductance, n, lo, hi))
        self.g_raw = nn.Parameter(_initial_raw(self.config.init_shunt_conductance, n, lo, hi))
        self.c_raw = nn.Parameter(_initial_raw(self.config.init_shunt_capacitance, n, lo, hi))
        self.training_history: list[dict[str, float]] = []

    def _param(self, raw: torch.Tensor) -> torch.Tensor:
        return _positive_range_from_raw(
            raw,
            self.config.parameter_min,
            self.config.parameter_max,
        )

    def normalized_parameters(self) -> dict[str, torch.Tensor]:
        """Return positive normalized Cauer element arrays."""

        return {
            "r": self._param(self.r_raw),
            "l": self._param(self.l_raw),
            "g": self._param(self.g_raw),
            "c": self._param(self.c_raw),
        }

    def initialize_from_model(self, source: "CauerLadderURN") -> None:
        """Copy a shallower ladder and append nearly-open tail sections."""

        if self.freqs_hz.shape != source.freqs_hz.shape or not np.allclose(
            self.freqs_hz, source.freqs_hz
        ):
            raise ValueError("source model must use the same frequency grid")
        if self.config.n_sections < source.config.n_sections:
            raise ValueError("target ladder must be at least as deep as source")

        lo = self.config.parameter_min
        hi = self.config.parameter_max
        source_params = source.normalized_parameters()
        n_copy = source.config.n_sections
        with torch.no_grad():
            self.r_raw.fill_(_raw_for_positive_value(self.config.extra_series_impedance, lo, hi))
            self.l_raw.fill_(_raw_for_positive_value(self.config.extra_series_impedance, lo, hi))
            self.g_raw.fill_(_raw_for_positive_value(self.config.extra_shunt_admittance, lo, hi))
            self.c_raw.fill_(_raw_for_positive_value(self.config.extra_shunt_admittance, lo, hi))
            self.r_raw[:n_copy].copy_(
                torch.tensor(
                    [
                        _raw_for_positive_value(float(x), lo, hi)
                        for x in source_params["r"].detach().cpu().numpy()
                    ],
                    dtype=torch.float64,
                    device=self.r_raw.device,
                )
            )
            self.l_raw[:n_copy].copy_(
                torch.tensor(
                    [
                        _raw_for_positive_value(float(x), lo, hi)
                        for x in source_params["l"].detach().cpu().numpy()
                    ],
                    dtype=torch.float64,
                    device=self.l_raw.device,
                )
            )
            self.g_raw[:n_copy].copy_(
                torch.tensor(
                    [
                        _raw_for_positive_value(float(x), lo, hi)
                        for x in source_params["g"].detach().cpu().numpy()
                    ],
                    dtype=torch.float64,
                    device=self.g_raw.device,
                )
            )
            self.c_raw[:n_copy].copy_(
                torch.tensor(
                    [
                        _raw_for_positive_value(float(x), lo, hi)
                        for x in source_params["c"].detach().cpu().numpy()
                    ],
                    dtype=torch.float64,
                    device=self.c_raw.device,
                )
            )

    def initialize_from_peeling(self, z_data: np.ndarray | list[complex]) -> None:
        """Initialize sections by alternating Z-domain and Y-domain peeling."""

        z_arr = _as_1d_complex128(z_data, "z_data")
        if z_arr.shape != self.freqs_hz.shape:
            raise ValueError("z_data must match this model's frequency grid")

        lo = self.config.parameter_min
        hi = self.config.parameter_max
        x_imag = self.omega / self.omega_ref
        z_work = z_arr / self.z_ref
        r_values: list[float] = []
        l_values: list[float] = []
        g_values: list[float] = []
        c_values: list[float] = []

        for _section in range(self.config.n_sections):
            r, l = _fit_positive_rx(z_work, x_imag)
            series = r + 1j * l * x_imag
            scale, z_residual = _stable_subtract(z_work, series)
            r *= scale
            l *= scale
            y_work = 1.0 / (z_residual + _EPS)

            g, c = _fit_positive_rx(y_work, x_imag)
            shunt = g + 1j * c * x_imag
            scale, y_residual = _stable_subtract(y_work, shunt)
            g *= scale
            c *= scale

            r_values.append(float(np.clip(r, lo, hi)))
            l_values.append(float(np.clip(l, lo, hi)))
            g_values.append(float(np.clip(g, lo, hi)))
            c_values.append(float(np.clip(c, lo, hi)))
            z_work = 1.0 / (y_residual + _EPS)

        with torch.no_grad():
            device = self.r_raw.device
            self.r_raw.copy_(
                torch.tensor(
                    [_raw_for_positive_value(x, lo, hi) for x in r_values],
                    dtype=torch.float64,
                    device=device,
                )
            )
            self.l_raw.copy_(
                torch.tensor(
                    [_raw_for_positive_value(x, lo, hi) for x in l_values],
                    dtype=torch.float64,
                    device=device,
                )
            )
            self.g_raw.copy_(
                torch.tensor(
                    [_raw_for_positive_value(x, lo, hi) for x in g_values],
                    dtype=torch.float64,
                    device=device,
                )
            )
            self.c_raw.copy_(
                torch.tensor(
                    [_raw_for_positive_value(x, lo, hi) for x in c_values],
                    dtype=torch.float64,
                    device=device,
                )
            )

    def _set_normalized_parameters(self, params: np.ndarray) -> None:
        params = np.asarray(params, dtype=np.float64).reshape(-1)
        n = self.config.n_sections
        if params.size != 4 * n:
            raise ValueError(f"expected {4 * n} Cauer parameters, got {params.size}")
        lo = self.config.parameter_min
        hi = self.config.parameter_max
        params = np.clip(params, lo, hi)
        with torch.no_grad():
            device = self.r_raw.device
            chunks = [params[i * n : (i + 1) * n] for i in range(4)]
            for raw, values in zip(
                (self.r_raw, self.l_raw, self.g_raw, self.c_raw),
                chunks,
                strict=True,
            ):
                raw.copy_(
                    torch.tensor(
                        [_raw_for_positive_value(float(x), lo, hi) for x in values],
                        dtype=torch.float64,
                        device=device,
                    )
                )

    def initialize_from_rational_fit(
        self,
        z_data: np.ndarray | list[complex],
        *,
        order: int | None = None,
        max_nfev: int | None = None,
    ) -> CauerRationalFit:
        """Initialize the positive Cauer ladder from a pole/zero rational teacher."""

        from scipy.optimize import least_squares

        z_arr = _as_1d_complex128(z_data, "z_data")
        if z_arr.shape != self.freqs_hz.shape:
            raise ValueError("z_data must match this model's frequency grid")

        rational = fit_rational_pole_zero(
            self.freqs_hz,
            z_arr,
            order=order or self.config.rational_order,
            omega_ref=self.omega_ref,
            z_ref=self.z_ref,
            sk_iterations=self.config.rational_sk_iterations,
            regularization=self.config.rational_regularization,
        )
        z_teacher = rational.predict(self.freqs_hz) / self.z_ref
        x = 1j * self.omega / self.omega_ref
        n = self.config.n_sections
        lo = self.config.parameter_min
        hi = self.config.parameter_max
        current = self.normalized_parameters()
        initial = np.concatenate(
            [
                current["r"].detach().cpu().numpy(),
                current["l"].detach().cpu().numpy(),
                current["g"].detach().cpu().numpy(),
                current["c"].detach().cpu().numpy(),
            ]
        )
        theta0 = np.log(np.clip(initial, lo, hi))
        lower = np.full(4 * n, np.log(lo), dtype=np.float64)
        upper = np.full(4 * n, np.log(hi), dtype=np.float64)

        def residual(theta: np.ndarray) -> np.ndarray:
            params = np.exp(theta)
            pred = _cauer_numpy_from_positive_params(params, x, n)
            return _relative_complex_residual(pred, z_teacher)

        result = least_squares(
            residual,
            theta0,
            bounds=(lower, upper),
            max_nfev=int(max_nfev or self.config.rational_cauer_max_nfev),
            x_scale="jac",
        )
        self._set_normalized_parameters(np.exp(result.x))
        self.initialization_info = {
            "rational_order": float(order or self.config.rational_order),
            "rational_cost": float(result.cost),
            "rational_nfev": float(result.nfev),
            "rational_success": float(bool(result.success)),
        }
        return rational

    def fit_to_response_least_squares(
        self,
        z_data: np.ndarray | list[complex],
        *,
        max_nfev: int | None = None,
        z_weight: float | None = None,
        y_weight: float | None = None,
    ) -> None:
        """Polish positive Cauer parameters against direct Z/Y residuals."""

        from scipy.optimize import least_squares

        z_arr = _as_1d_complex128(z_data, "z_data")
        if z_arr.shape != self.freqs_hz.shape:
            raise ValueError("z_data must match this model's frequency grid")

        n = self.config.n_sections
        lo = self.config.parameter_min
        hi = self.config.parameter_max
        x = 1j * self.omega / self.omega_ref
        z_target = z_arr / self.z_ref
        y_target = 1.0 / (z_target + _EPS)
        params = self.normalized_parameters()
        initial = np.concatenate(
            [
                params["r"].detach().cpu().numpy(),
                params["l"].detach().cpu().numpy(),
                params["g"].detach().cpu().numpy(),
                params["c"].detach().cpu().numpy(),
            ]
        )
        theta0 = np.log(np.clip(initial, lo, hi))
        lower = np.full(4 * n, np.log(lo), dtype=np.float64)
        upper = np.full(4 * n, np.log(hi), dtype=np.float64)
        zw = np.sqrt(self.config.least_squares_z_weight if z_weight is None else z_weight)
        yw = np.sqrt(self.config.least_squares_y_weight if y_weight is None else y_weight)

        def residual(theta: np.ndarray) -> np.ndarray:
            pred = _cauer_numpy_from_positive_params(np.exp(theta), x, n)
            y_pred = 1.0 / (pred + _EPS)
            return np.concatenate(
                [
                    zw * _relative_complex_residual(pred, z_target),
                    yw * _relative_complex_residual(y_pred, y_target),
                ]
            )

        result = least_squares(
            residual,
            theta0,
            bounds=(lower, upper),
            max_nfev=int(max_nfev or self.config.least_squares_max_nfev),
            x_scale="jac",
        )
        self._set_normalized_parameters(np.exp(result.x))
        self.least_squares_info = {
            "least_squares_cost": float(result.cost),
            "least_squares_nfev": float(result.nfev),
            "least_squares_success": float(bool(result.success)),
        }

    def normalized_impedance(self, omega: torch.Tensor) -> torch.Tensor:
        """Evaluate ``Z / z_ref`` by a Cauer continued fraction."""

        omega = omega.to(dtype=torch.float64).reshape(-1)
        x = 1j * omega / self.omega_ref
        params = self.normalized_parameters()
        r = params["r"]
        l = params["l"]
        g = params["g"]
        c = params["c"]
        z_next: torch.Tensor | None = None
        for i in range(self.config.n_sections - 1, -1, -1):
            y_tail = 0.0 if z_next is None else 1.0 / (z_next + _EPS)
            y_shunt = g[i] + c[i] * x + y_tail
            z_next = r[i] + l[i] * x + 1.0 / (y_shunt + _EPS)
        if z_next is None:
            raise RuntimeError("empty Cauer ladder")
        return z_next.to(torch.complex128)

    def forward(self, omega: torch.Tensor) -> torch.Tensor:
        """Evaluate the fitted impedance ``Z(omega)``."""

        z_ref = torch.as_tensor(self.z_ref, dtype=torch.float64, device=omega.device)
        return z_ref * self.normalized_impedance(omega)

    def _regularization(self) -> torch.Tensor:
        params = self.normalized_parameters()
        logs = [torch.log(v) for v in params.values()]
        reg = torch.mean(torch.cat([x * x for x in logs]))
        if self.config.log_smoothness_weight <= 0.0:
            return self.config.regularization_weight * reg
        smooth_terms = []
        for x in logs:
            if x.numel() > 1:
                smooth_terms.append(torch.mean(torch.diff(x) ** 2))
        smooth = (
            torch.sum(torch.stack(smooth_terms))
            if smooth_terms
            else torch.zeros((), dtype=torch.float64, device=self.r_raw.device)
        )
        return self.config.regularization_weight * reg + self.config.log_smoothness_weight * smooth

    def loss(self, omega: torch.Tensor, z_target: torch.Tensor) -> tuple[torch.Tensor, dict[str, float]]:
        z_pred = self(omega)
        residual = scattering_transform(z_pred, self.z0) - scattering_transform(z_target, self.z0)
        loss_fit = complex_smooth_l1(residual, self.config.huber_delta)
        loss_reg = self._regularization()
        loss = loss_fit + loss_reg
        return loss, {
            "loss": float(loss.detach().cpu()),
            "loss_fit": float(loss_fit.detach().cpu()),
            "loss_reg": float(loss_reg.detach().cpu()),
        }

    def zy_losses(
        self, omega: torch.Tensor, z_target: torch.Tensor
    ) -> tuple[torch.Tensor, torch.Tensor, torch.Tensor, dict[str, float]]:
        """Return direct impedance/admittance-domain losses.

        This avoids an S-parameter view.  Series-side parameters are naturally
        updated against the normalized impedance residual, while shunt-side
        parameters are naturally updated against the normalized admittance
        residual.
        """

        z_pred = self(omega)
        y_pred = 1.0 / (z_pred + _EPS)
        y_target = 1.0 / (z_target + _EPS)
        z_floor = torch.median(torch.abs(z_target)).clamp_min(1.0e-24) * 1.0e-6
        y_floor = torch.median(torch.abs(y_target)).clamp_min(1.0e-24) * 1.0e-6
        z_scale = torch.abs(z_target).clamp_min(z_floor)
        y_scale = torch.abs(y_target).clamp_min(y_floor)
        z_loss = complex_smooth_l1((z_pred - z_target) / z_scale, self.config.huber_delta)
        y_loss = complex_smooth_l1((y_pred - y_target) / y_scale, self.config.huber_delta)
        reg = self._regularization()
        total = z_loss + y_loss + reg
        return total, z_loss + reg, y_loss + reg, {
            "loss": float(total.detach().cpu()),
            "loss_z": float(z_loss.detach().cpu()),
            "loss_y": float(y_loss.detach().cpu()),
            "loss_reg": float(reg.detach().cpu()),
        }

    def predict(self, freqs_hz: np.ndarray | list[float]) -> np.ndarray:
        freqs = _as_1d_float64(freqs_hz, "freqs_hz")
        omega = torch.tensor(2.0 * np.pi * freqs, dtype=torch.float64, device=self.r_raw.device)
        with torch.no_grad():
            return self(omega).detach().cpu().numpy()

    def parameter_summary(self) -> list[dict[str, float | int]]:
        """Return normalized and physical section parameters."""

        params = {
            key: value.detach().cpu().numpy()
            for key, value in self.normalized_parameters().items()
        }
        out: list[dict[str, float | int]] = []
        for i in range(self.config.n_sections):
            r = float(params["r"][i])
            l = float(params["l"][i])
            g = float(params["g"][i])
            c = float(params["c"][i])
            out.append(
                {
                    "section": i,
                    "r_norm": r,
                    "l_norm": l,
                    "g_norm": g,
                    "c_norm": c,
                    "R_ohm": self.z_ref * r,
                    "L_h": self.z_ref * l / self.omega_ref,
                    "G_siemens": g / self.z_ref,
                    "C_f": c / (self.z_ref * self.omega_ref),
                }
            )
        return out


def train_cauer_ladder_urn(
    freqs_hz: np.ndarray | list[float],
    z_data: np.ndarray | list[complex],
    config: CauerLadderURNConfig | None = None,
    *,
    warm_start_model: CauerLadderURN | None = None,
    keep_best_epoch: bool = True,
    best_check_interval: int = 50,
    verbose: bool = True,
) -> CauerLadderURN:
    """Train one fixed-depth Cauer ladder with Adam."""

    freqs = _as_1d_float64(freqs_hz, "freqs_hz")
    z_arr = _as_1d_complex128(z_data, "z_data")
    if freqs.shape != z_arr.shape:
        raise ValueError("freqs_hz and z_data must have the same length")
    cfg = config or CauerLadderURNConfig()
    omega = torch.tensor(2.0 * np.pi * freqs, dtype=torch.float64)
    z_target = torch.tensor(z_arr, dtype=torch.complex128)

    best_model: CauerLadderURN | None = None
    best_loss = float("inf")
    best_history: list[dict[str, float]] = []
    check_interval = max(int(best_check_interval), 1)
    for restart_index in range(cfg.n_restarts):
        torch.manual_seed(int(cfg.seed) + 7919 * restart_index)
        model = CauerLadderURN(freqs, cfg, z_data=z_arr)
        if warm_start_model is not None:
            model.initialize_from_model(warm_start_model)
        elif cfg.use_rational_initialization:
            model.initialize_from_rational_fit(z_arr)
            if cfg.use_least_squares_polish:
                model.fit_to_response_least_squares(z_arr)
        elif cfg.use_peeling_initialization:
            model.initialize_from_peeling(z_arr)
        elif restart_index:
            with torch.no_grad():
                for param in (model.r_raw, model.l_raw, model.g_raw, model.c_raw):
                    param.add_(0.15 * torch.randn_like(param))

        optimizer = optim.Adam(model.parameters(), lr=cfg.lr)
        history: list[dict[str, float]] = []
        best_epoch_loss = float("inf")
        best_epoch_state: dict[str, torch.Tensor] | None = None

        def remember_if_best(epoch_number: int) -> None:
            nonlocal best_epoch_loss, best_epoch_state
            current_loss, current_parts = model.loss(omega, z_target)
            current_value = float(current_loss.detach().cpu())
            if current_value < best_epoch_loss:
                best_epoch_loss = current_value
                best_epoch_state = {
                    key: value.detach().clone()
                    for key, value in model.state_dict().items()
                }
            if epoch_number == 0 or epoch_number == cfg.n_epochs:
                parts = dict(current_parts)
                parts["epoch"] = float(epoch_number)
                history.append(parts)

        remember_if_best(0)
        for epoch in range(cfg.n_epochs):
            optimizer.zero_grad()
            loss, _parts = model.loss(omega, z_target)
            loss.backward()
            torch.nn.utils.clip_grad_norm_(model.parameters(), max_norm=10.0)
            optimizer.step()
            if (
                (keep_best_epoch and (epoch + 1) % check_interval == 0)
                or epoch == cfg.n_epochs - 1
            ):
                remember_if_best(epoch + 1)

        if keep_best_epoch and best_epoch_state is not None:
            model.load_state_dict(best_epoch_state)

        final_loss, parts = model.loss(omega, z_target)
        final_value = float(final_loss.detach().cpu())
        if verbose:
            print(
                f"  Cauer restart {restart_index + 1}/{cfg.n_restarts}: "
                f"loss={final_value:.6e}, fit={parts['loss_fit']:.6e}, "
                f"sections={cfg.n_sections}"
            )
        if final_value < best_loss:
            best_loss = final_value
            best_model = CauerLadderURN(freqs, cfg, z_data=z_arr)
            if warm_start_model is not None:
                best_model.initialize_from_model(warm_start_model)
            best_model.load_state_dict(model.state_dict())
            best_history = history

    if best_model is None:
        raise RuntimeError("Cauer ladder training did not produce a model")
    best_model.training_history = best_history
    return best_model


def train_cauer_ladder_alternating(
    freqs_hz: np.ndarray | list[float],
    z_data: np.ndarray | list[complex],
    config: CauerLadderURNConfig | None = None,
    *,
    warm_start_model: CauerLadderURN | None = None,
    active_sections: list[int] | range | None = None,
    best_check_interval: int = 1,
    verbose: bool = True,
) -> CauerLadderURN:
    """Train a Cauer ladder by alternating Z-side and Y-side updates.

    The series elements ``R,L`` are updated against a direct impedance residual.
    The shunt elements ``G,C`` are updated against a direct admittance residual.
    After each block, the combined Z+Y objective is checked; unstable blocks are
    rolled back and retried later with a reduced learning rate.
    """

    freqs = _as_1d_float64(freqs_hz, "freqs_hz")
    z_arr = _as_1d_complex128(z_data, "z_data")
    if freqs.shape != z_arr.shape:
        raise ValueError("freqs_hz and z_data must have the same length")
    cfg = config or CauerLadderURNConfig()
    omega = torch.tensor(2.0 * np.pi * freqs, dtype=torch.float64)
    z_target = torch.tensor(z_arr, dtype=torch.complex128)
    check_interval = max(int(best_check_interval), 1)
    if active_sections is None:
        active_section_list = list(range(cfg.n_sections))
    else:
        active_section_list = sorted({int(i) for i in active_sections})
    if not active_section_list:
        raise ValueError("active_sections must not be empty")
    if min(active_section_list) < 0 or max(active_section_list) >= cfg.n_sections:
        raise ValueError("active_sections contains an out-of-range section index")

    best_model: CauerLadderURN | None = None
    best_loss = float("inf")
    best_history: list[dict[str, float]] = []
    for restart_index in range(cfg.n_restarts):
        torch.manual_seed(int(cfg.seed) + 7919 * restart_index)
        model = CauerLadderURN(freqs, cfg, z_data=z_arr)
        if warm_start_model is not None:
            model.initialize_from_model(warm_start_model)
        elif cfg.use_rational_initialization:
            model.initialize_from_rational_fit(z_arr)
            if cfg.use_least_squares_polish:
                model.fit_to_response_least_squares(z_arr)
        elif cfg.use_peeling_initialization:
            model.initialize_from_peeling(z_arr)
        elif restart_index:
            with torch.no_grad():
                for param in (model.r_raw, model.l_raw, model.g_raw, model.c_raw):
                    param.add_(0.15 * torch.randn_like(param))

        lr_series = cfg.lr
        lr_shunt = cfg.lr
        history: list[dict[str, float]] = []
        best_epoch_state: dict[str, torch.Tensor] | None = None

        def objective_value() -> tuple[float, dict[str, float]]:
            total, _z_loss, _y_loss, parts = model.zy_losses(omega, z_target)
            return float(total.detach().cpu()), parts

        def snapshot() -> dict[str, torch.Tensor]:
            return {key: value.detach().clone() for key, value in model.state_dict().items()}

        def train_block(
            params: list[nn.Parameter],
            target: str,
            lr: float,
        ) -> tuple[float, bool]:
            before_value, _before_parts = objective_value()
            before_state = snapshot()
            optimizer = optim.Adam(params, lr=lr)
            for _ in range(cfg.alternating_block_epochs):
                optimizer.zero_grad()
                _total, z_loss, y_loss, _parts = model.zy_losses(omega, z_target)
                loss = z_loss if target == "z" else y_loss
                loss.backward()
                section_mask = torch.zeros(
                    cfg.n_sections,
                    dtype=torch.float64,
                    device=params[0].device,
                )
                section_mask[active_section_list] = 1.0
                for param in params:
                    if param.grad is not None:
                        param.grad.mul_(section_mask)
                torch.nn.utils.clip_grad_norm_(params, max_norm=10.0)
                optimizer.step()
            after_value, _after_parts = objective_value()
            if (
                not np.isfinite(after_value)
                or after_value > before_value * cfg.alternating_reject_growth
            ):
                model.load_state_dict(before_state)
                return before_value, False
            return after_value, True

        current_value, current_parts = objective_value()
        best_epoch_value = current_value
        best_epoch_state = snapshot()
        current_parts = dict(current_parts)
        current_parts["cycle"] = 0.0
        current_parts["accepted_series"] = 0.0
        current_parts["accepted_shunt"] = 0.0
        history.append(current_parts)

        series_failures = 0
        shunt_failures = 0
        for cycle in range(cfg.alternating_cycles):
            accepted_series = False
            accepted_shunt = False
            if lr_series >= cfg.alternating_min_lr:
                current_value, accepted_series = train_block(
                    [model.r_raw, model.l_raw],
                    "z",
                    lr_series,
                )
                if not accepted_series:
                    lr_series *= cfg.alternating_lr_decay_on_reject
                    series_failures += 1
            if lr_shunt >= cfg.alternating_min_lr:
                current_value, accepted_shunt = train_block(
                    [model.g_raw, model.c_raw],
                    "y",
                    lr_shunt,
                )
                if not accepted_shunt:
                    lr_shunt *= cfg.alternating_lr_decay_on_reject
                    shunt_failures += 1

            current_value, current_parts = objective_value()
            if current_value < best_epoch_value:
                best_epoch_value = current_value
                best_epoch_state = snapshot()
            if (cycle + 1) % check_interval == 0 or cycle == cfg.alternating_cycles - 1:
                parts = dict(current_parts)
                parts["cycle"] = float(cycle + 1)
                parts["accepted_series"] = float(accepted_series)
                parts["accepted_shunt"] = float(accepted_shunt)
                parts["lr_series"] = float(lr_series)
                parts["lr_shunt"] = float(lr_shunt)
                parts["series_failures"] = float(series_failures)
                parts["shunt_failures"] = float(shunt_failures)
                history.append(parts)

        if best_epoch_state is not None:
            model.load_state_dict(best_epoch_state)
        final_value, parts = objective_value()
        if verbose:
            print(
                f"  Cauer alternating restart {restart_index + 1}/{cfg.n_restarts}: "
                f"loss={final_value:.6e}, z={parts['loss_z']:.6e}, "
                f"y={parts['loss_y']:.6e}, sections={cfg.n_sections}"
            )
        if final_value < best_loss:
            best_loss = final_value
            best_model = CauerLadderURN(freqs, cfg, z_data=z_arr)
            if warm_start_model is not None:
                best_model.initialize_from_model(warm_start_model)
            best_model.load_state_dict(model.state_dict())
            best_history = history

    if best_model is None:
        raise RuntimeError("alternating Cauer ladder training did not produce a model")
    best_model.training_history = best_history
    return best_model


def train_cauer_ladder_tail_then_polish(
    freqs_hz: np.ndarray | list[float],
    z_data: np.ndarray | list[complex],
    config: CauerLadderURNConfig | None = None,
    *,
    frozen_outer_sections: int | None = None,
    verbose: bool = True,
) -> CauerLadderURN:
    """Train inner Cauer sections with outer sections frozen, then polish all."""

    cfg = config or CauerLadderURNConfig.twenty_two_parameter_candidate()
    frozen = cfg.frozen_outer_sections if frozen_outer_sections is None else int(frozen_outer_sections)
    frozen = max(0, min(frozen, cfg.n_sections - 1))
    if frozen == 0:
        return train_cauer_ladder_alternating(freqs_hz, z_data, cfg, verbose=verbose)

    tail_cycles = cfg.tail_train_cycles or max(1, cfg.alternating_cycles // 2)
    polish_cycles = cfg.polish_train_cycles or cfg.alternating_cycles
    tail_cfg = replace(cfg, alternating_cycles=tail_cycles)
    polish_cfg = replace(cfg, alternating_cycles=polish_cycles)
    active_tail = range(frozen, cfg.n_sections)
    if verbose:
        print(
            f"Cauer tail stage: freeze outer {frozen} section(s), "
            f"train sections {frozen}..{cfg.n_sections - 1}"
        )
    tail_model = train_cauer_ladder_alternating(
        freqs_hz,
        z_data,
        tail_cfg,
        active_sections=active_tail,
        verbose=verbose,
    )
    if verbose:
        print("Cauer polish stage: unfreeze all sections")
    model = train_cauer_ladder_alternating(
        freqs_hz,
        z_data,
        polish_cfg,
        warm_start_model=tail_model,
        active_sections=None,
        verbose=verbose,
    )
    model.training_history = [
        {"stage": "tail", **item} for item in tail_model.training_history
    ] + [{"stage": "polish", **item} for item in model.training_history]
    return model


def train_cauer_ladder_progressive(
    freqs_hz: np.ndarray | list[float],
    z_data: np.ndarray | list[complex],
    config: CauerLadderURNConfig | None = None,
    *,
    method: str = "alternating",
    verbose: bool = True,
) -> CauerLadderURN:
    """Fit a Cauer ladder by growing the continued fraction one section at a time."""

    cfg = config or CauerLadderURNConfig.twenty_two_parameter_candidate()
    start = max(1, min(int(cfg.progressive_start_sections), cfg.n_sections))
    stage_epochs = cfg.progressive_stage_epochs or cfg.n_epochs
    model: CauerLadderURN | None = None
    stage_history: list[dict[str, float]] = []
    for stage_index, n_sections in enumerate(range(start, cfg.n_sections + 1)):
        stage_cfg = replace(
            cfg,
            n_sections=n_sections,
            n_epochs=stage_epochs,
            lr=cfg.lr * (cfg.progressive_lr_decay ** stage_index),
        )
        if verbose:
            print(
                f"Cauer stage {stage_index + 1}: "
                f"sections={n_sections}, parameters={stage_cfg.total_parameters}"
            )
        if method == "adam":
            model = train_cauer_ladder_urn(
                freqs_hz,
                z_data,
                stage_cfg,
                warm_start_model=model,
                verbose=verbose,
            )
        elif method == "alternating":
            model = train_cauer_ladder_alternating(
                freqs_hz,
                z_data,
                stage_cfg,
                warm_start_model=model,
                verbose=verbose,
            )
        else:
            raise ValueError("method must be 'alternating' or 'adam'")
        if model.training_history:
            stage_history.append(
                {
                    "stage": float(stage_index + 1),
                    "sections": float(n_sections),
                    **model.training_history[-1],
                }
            )
    if model is None:
        raise RuntimeError("progressive Cauer ladder training did not produce a model")
    model.training_history = stage_history
    return model


__all__ = [
    "CauerLadderURN",
    "CauerLadderURNConfig",
    "CauerRationalFit",
    "fit_rational_pole_zero",
    "train_cauer_ladder_alternating",
    "train_cauer_ladder_progressive",
    "train_cauer_ladder_tail_then_polish",
    "train_cauer_ladder_urn",
]
