"""Y-domain Universal Relaxation Network variant.

This module implements the admittance-domain formulation used in the
SA/RM-2026 URN manuscript draft: physical basis functions are summed in
admittance space, fitted through an S-domain loss, and pruned by output
ablation.  It is intentionally additive to the older ``UniversalRelaxationNetwork``
API, whose series/parallel Z-domain model remains the default production path.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any

import numpy as np
import torch
from torch import nn, optim


_EPS = 1.0e-30


@dataclass
class YAdmittanceURNConfig:
    """Configuration for the SA/RM-2026 Y-base URN.

    The defaults reproduce the 22-basis dictionary described in the research
    meeting draft:

    Debye 5, magnetic Debye 5, Cole-Cole 3, magnetic Cole-Cole 3,
    inductive CPE 2, capacitive CPE 2, and series-RLC resonance 2.
    """

    n_debye: int = 5
    n_magnetic_debye: int = 5
    n_cole_cole: int = 3
    n_magnetic_cole_cole: int = 3
    n_inductive_cpe: int = 2
    n_capacitive_cpe: int = 2
    n_series_rlc: int = 2
    # Anti-resonance extensions (0 by default = the paper 22-basis dictionary).
    # parallel_rlc: R||L||C, Y*R = 1 + jQ(x - 1/x)  (damped anti-resonance)
    # coil_antiresonance: C||(R+L), Y*R = jx/Q + 1/(1+jQx)  (self-resonant coil)
    n_parallel_rlc: int = 0
    n_coil_antiresonance: int = 0

    sparsity_weight: float = 1.0e-4
    lr: float = 1.0e-3
    n_epochs: int = 5000
    n_restarts: int = 1
    seed: int = 7

    huber_delta: float = 1.0e-2
    active_threshold: float = 1.0e-3
    gate_init: float = 0.05
    alpha_min: float = 0.001
    alpha_max: float = 0.95
    beta_min: float = 0.1
    beta_max: float = 0.95
    q_min: float = 0.05
    q_max: float = 20.0
    tau_min: float | None = None
    tau_max: float | None = None
    omega_ref: float | None = None
    z0: float | None = None
    y_ref: float | None = None

    @property
    def total_basis_functions(self) -> int:
        return (
            self.n_debye
            + self.n_magnetic_debye
            + self.n_cole_cole
            + self.n_magnetic_cole_cole
            + self.n_inductive_cpe
            + self.n_capacitive_cpe
            + self.n_series_rlc
            + self.n_parallel_rlc
            + self.n_coil_antiresonance
        )

    @classmethod
    def paper_22_basis(cls, **overrides: Any) -> "YAdmittanceURNConfig":
        """Return the research-meeting 22-basis attention-free configuration."""

        values = {
            "sparsity_weight": 1.0e-4,
            "lr": 5.0e-3,
            "n_epochs": 50_000,
            "n_restarts": 1,
            "gate_init": 1.0,
            "active_threshold": 1.0e-3,
        }
        values.update(overrides)
        return cls(**values)


@dataclass(frozen=True)
class YAdmittanceURNActiveBasis:
    """Output-ablation summary for one active Y-base URN basis."""

    basis_index: int
    basis_type: str
    local_index: int
    importance: float
    gate: float
    parameters: dict[str, float] = field(default_factory=dict)


def _as_1d_float64(values: np.ndarray | list[float], name: str) -> np.ndarray:
    arr = np.asarray(values, dtype=np.float64).reshape(-1)
    if arr.size == 0:
        raise ValueError(f"{name} must not be empty")
    if not np.all(np.isfinite(arr)):
        raise ValueError(f"{name} contains non-finite values")
    return arr


def _as_1d_complex128(values: np.ndarray | list[complex], name: str) -> np.ndarray:
    arr = np.asarray(values, dtype=np.complex128).reshape(-1)
    if arr.size == 0:
        raise ValueError(f"{name} must not be empty")
    if not np.all(np.isfinite(arr.real)) or not np.all(np.isfinite(arr.imag)):
        raise ValueError(f"{name} contains non-finite values")
    return arr


def _inv_sigmoid(x: np.ndarray | float) -> np.ndarray:
    x_arr = np.asarray(x, dtype=np.float64)
    x_arr = np.clip(x_arr, 1.0e-6, 1.0 - 1.0e-6)
    return np.log(x_arr / (1.0 - x_arr))


def _bounded_from_raw(raw: torch.Tensor, lo: float, hi: float) -> torch.Tensor:
    return lo + (hi - lo) * torch.sigmoid(raw)


def _positive_range_from_raw(raw: torch.Tensor, lo: float, hi: float) -> torch.Tensor:
    log_lo = torch.as_tensor(np.log(lo), dtype=torch.float64, device=raw.device)
    log_hi = torch.as_tensor(np.log(hi), dtype=torch.float64, device=raw.device)
    return torch.exp(log_lo + (log_hi - log_lo) * torch.sigmoid(raw))


def _initial_raw_log_grid(n: int, lo: float, hi: float) -> torch.Tensor:
    if n <= 0:
        return torch.empty(0, dtype=torch.float64)
    if n == 1:
        fractions = np.array([0.5], dtype=np.float64)
    else:
        # Stay away from the exact bounds so optimization can move both ways.
        fractions = np.linspace(0.1, 0.9, n, dtype=np.float64)
    return torch.tensor(_inv_sigmoid(fractions), dtype=torch.float64)


def _initial_raw_constant(value: float, lo: float, hi: float, n: int) -> torch.Tensor:
    if n <= 0:
        return torch.empty(0, dtype=torch.float64)
    frac = (value - lo) / (hi - lo)
    return torch.full((n,), float(_inv_sigmoid(frac)), dtype=torch.float64)


def _raw_for_positive_value(value: float, lo: float, hi: float) -> float:
    value = float(np.clip(value, lo * (1.0 + 1.0e-9), hi * (1.0 - 1.0e-9)))
    frac = (np.log(value) - np.log(lo)) / (np.log(hi) - np.log(lo))
    return float(_inv_sigmoid(frac))


def _raw_for_bounded_value(value: float, lo: float, hi: float) -> float:
    frac = (float(value) - lo) / (hi - lo)
    return float(_inv_sigmoid(frac))


def _raw_for_gate(value: float) -> float:
    gate = max(float(value), 1.0e-9)
    return float(np.log(np.expm1(gate)))


def scattering_transform(z: torch.Tensor, z0: float) -> torch.Tensor:
    """Map impedance to a bounded S-like response using a real reference z0."""

    z0_t = torch.as_tensor(float(z0), dtype=torch.float64, device=z.device)
    return (z - z0_t) / (z + z0_t + _EPS)


def s_domain_rmse(
    z_fit: np.ndarray | list[complex],
    z_target: np.ndarray | list[complex],
    *,
    z0: float | None = None,
) -> float:
    """Return the S-domain RMSE used by the SA/RM-2026 URN manuscript."""

    fit = _as_1d_complex128(z_fit, "z_fit")
    target = _as_1d_complex128(z_target, "z_target")
    if fit.shape != target.shape:
        raise ValueError("z_fit and z_target must have the same length")
    z0_value = float(np.median(np.abs(target))) if z0 is None else float(z0)
    z0_value = max(z0_value, _EPS)
    s_fit = (fit - z0_value) / (fit + z0_value + _EPS)
    s_target = (target - z0_value) / (target + z0_value + _EPS)
    return float(np.sqrt(np.mean(np.abs(s_fit - s_target) ** 2)))


def complex_smooth_l1(
    residual: torch.Tensor,
    beta: float,
    weight: torch.Tensor | None = None,
) -> torch.Tensor:
    """Smooth-L1 loss on real and imaginary residual components.

    ``weight`` optionally supplies nonnegative per-sample weights; the result
    is then the weighted mean ``sum(w * l) / sum(w)``.
    """

    def one_part(x: torch.Tensor) -> torch.Tensor:
        ax = torch.abs(x)
        beta_t = torch.as_tensor(float(beta), dtype=torch.float64, device=x.device)
        return torch.where(ax < beta_t, 0.5 * ax * ax / beta_t, ax - 0.5 * beta_t)

    per_point = one_part(residual.real) + one_part(residual.imag)
    if weight is None:
        return torch.mean(per_point)
    w = weight.to(dtype=per_point.dtype, device=per_point.device)
    return torch.sum(w * per_point) / torch.sum(w).clamp_min(_EPS)


_BASIS_COUNT_FIELDS = {
    "debye": "n_debye",
    "magnetic_debye": "n_magnetic_debye",
    "cole_cole": "n_cole_cole",
    "magnetic_cole_cole": "n_magnetic_cole_cole",
    "inductive_cpe": "n_inductive_cpe",
    "capacitive_cpe": "n_capacitive_cpe",
    "series_rlc": "n_series_rlc",
    "parallel_rlc": "n_parallel_rlc",
    "coil_antiresonance": "n_coil_antiresonance",
}


def active_basis_refit_config(
    active_bases: list[YAdmittanceURNActiveBasis],
    base_config: YAdmittanceURNConfig | None = None,
    **overrides: Any,
) -> YAdmittanceURNConfig:
    """Build an attention-free config containing only selected active bases."""

    active = list(active_bases)
    if not active:
        raise ValueError("active_bases must not be empty")

    counts = dict.fromkeys(_BASIS_COUNT_FIELDS.values(), 0)
    for basis in active:
        if basis.basis_type not in _BASIS_COUNT_FIELDS:
            raise ValueError(f"unsupported active basis type: {basis.basis_type!r}")
        counts[_BASIS_COUNT_FIELDS[basis.basis_type]] += 1

    base = base_config or YAdmittanceURNConfig.paper_22_basis()
    values = {
        field_name: getattr(base, field_name)
        for field_name in YAdmittanceURNConfig.__dataclass_fields__
    }
    values.update(counts)
    values.update(
        {
            "sparsity_weight": 0.0,
            "gate_init": 1.0,
        }
    )
    values.update(overrides)
    return YAdmittanceURNConfig(**values)


class YAdmittanceURN(nn.Module):
    """Admittance-domain URN with nonnegative frequency-independent gates."""

    def __init__(
        self,
        freqs_hz: np.ndarray | list[float],
        config: YAdmittanceURNConfig | None = None,
        *,
        z_data: np.ndarray | list[complex] | None = None,
    ) -> None:
        super().__init__()
        self.config = config or YAdmittanceURNConfig()
        if self.config.total_basis_functions <= 0:
            raise ValueError("at least one basis function is required")

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

        tau_min = self.config.tau_min or (0.1 / self.omega_max)
        tau_max = self.config.tau_max or (10.0 / self.omega_min)
        if tau_min <= 0.0 or tau_max <= tau_min:
            raise ValueError("tau_min/tau_max must be positive and ordered")
        self.tau_min = float(tau_min)
        self.tau_max = float(tau_max)

        if z_data is not None:
            z_arr = _as_1d_complex128(z_data, "z_data")
            if z_arr.shape != freqs.shape:
                raise ValueError("freqs_hz and z_data must have the same length")
            y_arr = 1.0 / z_arr
            if not np.all(np.isfinite(y_arr.real)) or not np.all(np.isfinite(y_arr.imag)):
                raise ValueError("z_data contains zero or singular impedance values")
            z0 = self.config.z0 or float(np.median(np.abs(z_arr)))
            y_ref = self.config.y_ref or float(np.mean(np.abs(y_arr)))
        else:
            z0 = self.config.z0 or 1.0
            y_ref = self.config.y_ref or 1.0
        self.z0 = float(max(z0, _EPS))
        self.y_ref = float(max(y_ref, _EPS))

        self.basis_labels = self._build_basis_labels()
        n_basis = len(self.basis_labels)
        gate_init = max(float(self.config.gate_init), 1.0e-9)
        self.gate_raw = nn.Parameter(
            torch.full((n_basis,), np.log(np.expm1(gate_init)), dtype=torch.float64)
        )

        self._init_parameters()

        self.training_history: list[dict[str, float]] = []

    def _build_basis_labels(self) -> list[tuple[str, int]]:
        labels: list[tuple[str, int]] = []
        for name, count in (
            ("debye", self.config.n_debye),
            ("magnetic_debye", self.config.n_magnetic_debye),
            ("cole_cole", self.config.n_cole_cole),
            ("magnetic_cole_cole", self.config.n_magnetic_cole_cole),
            ("inductive_cpe", self.config.n_inductive_cpe),
            ("capacitive_cpe", self.config.n_capacitive_cpe),
            ("series_rlc", self.config.n_series_rlc),
            ("parallel_rlc", self.config.n_parallel_rlc),
            ("coil_antiresonance", self.config.n_coil_antiresonance),
        ):
            labels.extend((name, i) for i in range(count))
        return labels

    def _init_parameters(self) -> None:
        c = self.config
        self.debye_tau_raw = nn.Parameter(_initial_raw_log_grid(c.n_debye, self.tau_min, self.tau_max))
        self.mag_debye_tau_raw = nn.Parameter(
            _initial_raw_log_grid(c.n_magnetic_debye, self.tau_min, self.tau_max)
        )
        self.cc_tau_raw = nn.Parameter(_initial_raw_log_grid(c.n_cole_cole, self.tau_min, self.tau_max))
        self.cc_alpha_raw = nn.Parameter(
            _initial_raw_constant(0.75, c.alpha_min, c.alpha_max, c.n_cole_cole)
        )
        self.mag_cc_tau_raw = nn.Parameter(
            _initial_raw_log_grid(c.n_magnetic_cole_cole, self.tau_min, self.tau_max)
        )
        self.mag_cc_alpha_raw = nn.Parameter(
            _initial_raw_constant(0.75, c.alpha_min, c.alpha_max, c.n_magnetic_cole_cole)
        )
        self.ind_cpe_beta_raw = nn.Parameter(
            _initial_raw_constant(0.9, c.beta_min, c.beta_max, c.n_inductive_cpe)
        )
        self.cap_cpe_beta_raw = nn.Parameter(
            _initial_raw_constant(0.9, c.beta_min, c.beta_max, c.n_capacitive_cpe)
        )
        self.rlc_omega0_raw = nn.Parameter(
            _initial_raw_log_grid(c.n_series_rlc, self.omega_min, self.omega_max)
        )
        self.rlc_q_raw = nn.Parameter(_initial_raw_constant(0.7, c.q_min, c.q_max, c.n_series_rlc))
        self.prlc_omega0_raw = nn.Parameter(
            _initial_raw_log_grid(c.n_parallel_rlc, self.omega_min, self.omega_max)
        )
        self.prlc_q_raw = nn.Parameter(
            _initial_raw_constant(0.7, c.q_min, c.q_max, c.n_parallel_rlc)
        )
        self.coil_omega0_raw = nn.Parameter(
            _initial_raw_log_grid(c.n_coil_antiresonance, self.omega_min, self.omega_max)
        )
        self.coil_q_raw = nn.Parameter(
            _initial_raw_constant(2.0, c.q_min, c.q_max, c.n_coil_antiresonance)
        )

    def initialize_from_active_bases(
        self, active_bases: list[YAdmittanceURNActiveBasis]
    ) -> None:
        """Initialize this attention-free model from selected active bases."""

        active = list(active_bases)
        if len(active) != len(self.basis_labels):
            raise ValueError(
                "active_bases count must match this model's basis count; "
                f"got {len(active)} for {len(self.basis_labels)} bases"
            )

        by_type: dict[str, list[YAdmittanceURNActiveBasis]] = {
            name: [] for name in _BASIS_COUNT_FIELDS
        }
        for basis in active:
            if basis.basis_type not in by_type:
                raise ValueError(f"unsupported active basis type: {basis.basis_type!r}")
            by_type[basis.basis_type].append(basis)
        for basis_type, config_field in _BASIS_COUNT_FIELDS.items():
            expected = getattr(self.config, config_field)
            if len(by_type[basis_type]) != expected:
                raise ValueError(
                    f"{basis_type!r} active basis count {len(by_type[basis_type])} "
                    f"does not match config count {expected}"
                )

        def require(basis: YAdmittanceURNActiveBasis, key: str) -> float:
            if key not in basis.parameters:
                raise ValueError(f"{basis.basis_type} active basis is missing parameter {key!r}")
            return float(basis.parameters[key])

        def tau_raw(items: list[YAdmittanceURNActiveBasis]) -> torch.Tensor:
            return torch.tensor(
                [
                    _raw_for_positive_value(require(item, "tau"), self.tau_min, self.tau_max)
                    for item in items
                ],
                dtype=torch.float64,
            )

        def bounded_raw(
            items: list[YAdmittanceURNActiveBasis], key: str, lo: float, hi: float
        ) -> torch.Tensor:
            return torch.tensor(
                [_raw_for_bounded_value(require(item, key), lo, hi) for item in items],
                dtype=torch.float64,
            )

        with torch.no_grad():
            ordered = [by_type[basis_type][local_index] for basis_type, local_index in self.basis_labels]
            self.gate_raw.copy_(
                torch.tensor([_raw_for_gate(item.gate) for item in ordered], dtype=torch.float64)
            )
            self.debye_tau_raw.copy_(tau_raw(by_type["debye"]))
            self.mag_debye_tau_raw.copy_(tau_raw(by_type["magnetic_debye"]))
            self.cc_tau_raw.copy_(tau_raw(by_type["cole_cole"]))
            self.cc_alpha_raw.copy_(
                bounded_raw(
                    by_type["cole_cole"],
                    "alpha",
                    self.config.alpha_min,
                    self.config.alpha_max,
                )
            )
            self.mag_cc_tau_raw.copy_(tau_raw(by_type["magnetic_cole_cole"]))
            self.mag_cc_alpha_raw.copy_(
                bounded_raw(
                    by_type["magnetic_cole_cole"],
                    "alpha",
                    self.config.alpha_min,
                    self.config.alpha_max,
                )
            )
            self.ind_cpe_beta_raw.copy_(
                bounded_raw(
                    by_type["inductive_cpe"],
                    "beta",
                    self.config.beta_min,
                    self.config.beta_max,
                )
            )
            self.cap_cpe_beta_raw.copy_(
                bounded_raw(
                    by_type["capacitive_cpe"],
                    "beta",
                    self.config.beta_min,
                    self.config.beta_max,
                )
            )
            def omega0_raw(basis_type: str) -> torch.Tensor:
                omega0_values = []
                for item in by_type[basis_type]:
                    if "omega0" in item.parameters:
                        omega0_values.append(float(item.parameters["omega0"]))
                    elif "f0" in item.parameters:
                        omega0_values.append(2.0 * np.pi * float(item.parameters["f0"]))
                    else:
                        raise ValueError(
                            f"{basis_type} active basis is missing 'omega0' or 'f0'"
                        )
                return torch.tensor(
                    [
                        _raw_for_positive_value(value, self.omega_min, self.omega_max)
                        for value in omega0_values
                    ],
                    dtype=torch.float64,
                )

            self.rlc_omega0_raw.copy_(omega0_raw("series_rlc"))
            self.rlc_q_raw.copy_(
                bounded_raw(by_type["series_rlc"], "q", self.config.q_min, self.config.q_max)
            )
            self.prlc_omega0_raw.copy_(omega0_raw("parallel_rlc"))
            self.prlc_q_raw.copy_(
                bounded_raw(by_type["parallel_rlc"], "q", self.config.q_min, self.config.q_max)
            )
            self.coil_omega0_raw.copy_(omega0_raw("coil_antiresonance"))
            self.coil_q_raw.copy_(
                bounded_raw(
                    by_type["coil_antiresonance"], "q", self.config.q_min, self.config.q_max
                )
            )

    def initialize_from_model(
        self,
        source: "YAdmittanceURN",
        *,
        extra_gate: float = 1.0e-8,
    ) -> None:
        """Embed a smaller fitted dictionary into this model.

        Matching basis labels keep their fitted gates and physical parameters.
        Bases that exist only in this model are initialized with a nearly-zero
        gate, so a larger dictionary starts from the smaller model's response.
        """

        if self.freqs_hz.shape != source.freqs_hz.shape or not np.allclose(
            self.freqs_hz, source.freqs_hz
        ):
            raise ValueError("source model must use the same frequency grid")

        source_label_to_index = {label: i for i, label in enumerate(source.basis_labels)}
        source_gates = source.gates().detach().cpu().numpy()
        with torch.no_grad():
            self.gate_raw.fill_(_raw_for_gate(extra_gate))
            for target_index, label in enumerate(self.basis_labels):
                source_index = source_label_to_index.get(label)
                if source_index is None:
                    continue
                self.gate_raw[target_index] = _raw_for_gate(source_gates[source_index])

            def copy_prefix(target: torch.Tensor, src: torch.Tensor) -> None:
                n = min(target.numel(), src.numel())
                if n:
                    target[:n].copy_(src[:n])

            copy_prefix(self.debye_tau_raw, source.debye_tau_raw)
            copy_prefix(self.mag_debye_tau_raw, source.mag_debye_tau_raw)
            copy_prefix(self.cc_tau_raw, source.cc_tau_raw)
            copy_prefix(self.cc_alpha_raw, source.cc_alpha_raw)
            copy_prefix(self.mag_cc_tau_raw, source.mag_cc_tau_raw)
            copy_prefix(self.mag_cc_alpha_raw, source.mag_cc_alpha_raw)
            copy_prefix(self.ind_cpe_beta_raw, source.ind_cpe_beta_raw)
            copy_prefix(self.cap_cpe_beta_raw, source.cap_cpe_beta_raw)
            copy_prefix(self.rlc_omega0_raw, source.rlc_omega0_raw)
            copy_prefix(self.rlc_q_raw, source.rlc_q_raw)
            copy_prefix(self.prlc_omega0_raw, source.prlc_omega0_raw)
            copy_prefix(self.prlc_q_raw, source.prlc_q_raw)
            copy_prefix(self.coil_omega0_raw, source.coil_omega0_raw)
            copy_prefix(self.coil_q_raw, source.coil_q_raw)

    def _tau(self, raw: torch.Tensor) -> torch.Tensor:
        return _positive_range_from_raw(raw, self.tau_min, self.tau_max)

    def _omega0(self, raw: torch.Tensor) -> torch.Tensor:
        return _positive_range_from_raw(raw, self.omega_min, self.omega_max)

    def gates(self) -> torch.Tensor:
        """Return nonnegative dimensionless gates."""

        return torch.nn.functional.softplus(self.gate_raw)

    def basis_matrix(self, omega: torch.Tensor, *, normalize: bool = True) -> torch.Tensor:
        """Return basis responses in admittance space with shape ``(N, M)``."""

        omega = omega.to(dtype=torch.float64).reshape(-1)
        jw = 1j * omega
        w_norm = omega / self.omega_ref
        jw_norm = 1j * w_norm
        cols: list[torch.Tensor] = []

        if self.config.n_debye:
            tau = self._tau(self.debye_tau_raw)
            x = jw[:, None] * tau[None, :]
            cols.append(1.0 / (1.0 + x))

        if self.config.n_magnetic_debye:
            tau = self._tau(self.mag_debye_tau_raw)
            x = jw[:, None] * tau[None, :]
            cols.append(x / (1.0 + x))

        if self.config.n_cole_cole:
            tau = self._tau(self.cc_tau_raw)
            alpha = _bounded_from_raw(
                self.cc_alpha_raw, self.config.alpha_min, self.config.alpha_max
            )
            x = (jw[:, None] * tau[None, :]) ** alpha[None, :]
            cols.append(1.0 / (1.0 + x))

        if self.config.n_magnetic_cole_cole:
            tau = self._tau(self.mag_cc_tau_raw)
            alpha = _bounded_from_raw(
                self.mag_cc_alpha_raw, self.config.alpha_min, self.config.alpha_max
            )
            x_linear = jw[:, None] * tau[None, :]
            x_frac = x_linear ** alpha[None, :]
            cols.append(x_linear / (1.0 + x_frac))

        if self.config.n_inductive_cpe:
            beta = _bounded_from_raw(
                self.ind_cpe_beta_raw, self.config.beta_min, self.config.beta_max
            )
            cols.append(jw_norm[:, None] ** (-beta[None, :]))

        if self.config.n_capacitive_cpe:
            beta = _bounded_from_raw(
                self.cap_cpe_beta_raw, self.config.beta_min, self.config.beta_max
            )
            cols.append(jw_norm[:, None] ** beta[None, :])

        if self.config.n_series_rlc:
            omega0 = self._omega0(self.rlc_omega0_raw)
            q = _bounded_from_raw(self.rlc_q_raw, self.config.q_min, self.config.q_max)
            ratio = omega[:, None] / omega0[None, :]
            cols.append(1.0 / (1.0 + 1j * q[None, :] * (ratio - 1.0 / ratio)))

        if self.config.n_parallel_rlc:
            # R||L||C: Y*R = 1 + jQ(x - 1/x); Re(Y)*R = 1 > 0 (damped
            # anti-resonance -- admittance valley / impedance peak at x=1).
            omega0 = self._omega0(self.prlc_omega0_raw)
            q = _bounded_from_raw(self.prlc_q_raw, self.config.q_min, self.config.q_max)
            ratio = omega[:, None] / omega0[None, :]
            cols.append(
                (1.0 + 1j * q[None, :] * (ratio - 1.0 / ratio)).to(torch.complex128)
            )

        if self.config.n_coil_antiresonance:
            # C||(R+L) self-resonant coil: Y*R = jx/Q + 1/(1 + jQx);
            # Re(Y)*R = 1/(1+Q^2 x^2) > 0 (damped anti-resonance near x=1).
            omega0 = self._omega0(self.coil_omega0_raw)
            q = _bounded_from_raw(self.coil_q_raw, self.config.q_min, self.config.q_max)
            ratio = omega[:, None] / omega0[None, :]
            cols.append(
                1j * ratio / q[None, :] + 1.0 / (1.0 + 1j * q[None, :] * ratio)
            )

        if not cols:
            raise RuntimeError("basis dictionary is empty")

        basis = torch.cat(cols, dim=1).to(torch.complex128)
        if normalize:
            rms = torch.sqrt(torch.mean(torch.abs(basis) ** 2, dim=0)).clamp_min(1.0e-24)
            basis = basis / rms[None, :]
        return basis

    def forward_admittance(self, omega: torch.Tensor, mask: torch.Tensor | None = None) -> torch.Tensor:
        """Evaluate the fitted admittance ``Y(omega)``."""

        basis = self.basis_matrix(omega, normalize=True)
        weights = self.gates()[None, :]
        if mask is not None:
            weights = weights * mask.to(dtype=torch.float64, device=weights.device)[None, :]
        return self.y_ref * torch.sum(weights.to(torch.complex128) * basis, dim=1)

    def forward(self, omega: torch.Tensor, mask: torch.Tensor | None = None) -> torch.Tensor:
        """Evaluate the fitted impedance ``Z(omega)=1/Y(omega)``."""

        y = self.forward_admittance(omega, mask=mask)
        y_floor = torch.as_tensor(self.y_ref * 1.0e-24, dtype=torch.float64, device=y.device)
        return 1.0 / (y + y_floor)

    def loss(self, omega: torch.Tensor, z_target: torch.Tensor) -> tuple[torch.Tensor, dict[str, float]]:
        z_pred = self(omega)
        residual = scattering_transform(z_pred, self.z0) - scattering_transform(z_target, self.z0)
        loss_fit = complex_smooth_l1(residual, self.config.huber_delta)
        loss_sparse = self.config.sparsity_weight * torch.mean(self.gates())
        loss = loss_fit + loss_sparse
        return loss, {
            "loss": float(loss.detach().cpu()),
            "loss_fit": float(loss_fit.detach().cpu()),
            "loss_sparse": float(loss_sparse.detach().cpu()),
            "mean_gate": float(torch.mean(self.gates()).detach().cpu()),
        }

    def predict(self, freqs_hz: np.ndarray | list[float]) -> np.ndarray:
        freqs = _as_1d_float64(freqs_hz, "freqs_hz")
        omega = torch.tensor(2.0 * np.pi * freqs, dtype=torch.float64)
        with torch.no_grad():
            return self(omega).detach().cpu().numpy()

    def output_ablation_importance(
        self, freqs_hz: np.ndarray | list[float] | None = None
    ) -> np.ndarray:
        """Return output-ablation importance ``I_i`` for every basis."""

        freqs = self.freqs_hz if freqs_hz is None else _as_1d_float64(freqs_hz, "freqs_hz")
        omega = torch.tensor(2.0 * np.pi * freqs, dtype=torch.float64)
        n_basis = len(self.basis_labels)
        with torch.no_grad():
            z_all = self(omega)
            denom = torch.sqrt(torch.mean(torch.abs(z_all) ** 2)).clamp_min(1.0e-24)
            out = []
            for i in range(n_basis):
                mask = torch.ones(n_basis, dtype=torch.float64)
                mask[i] = 0.0
                z_without = self(omega, mask=mask)
                importance = torch.sqrt(torch.mean(torch.abs(z_all - z_without) ** 2)) / denom
                out.append(float(importance.detach().cpu()))
        return np.asarray(out, dtype=np.float64)

    def parameter_summary(self) -> list[dict[str, float | int | str]]:
        """Return current physical parameters in basis order."""

        summaries: list[dict[str, float | int | str]] = []
        gates = self.gates().detach().cpu().numpy()
        param_iters = self._parameter_values_by_type()
        for basis_index, (basis_type, local_index) in enumerate(self.basis_labels):
            params = dict(param_iters[basis_type][local_index])
            summaries.append(
                {
                    "basis_index": basis_index,
                    "basis_type": basis_type,
                    "local_index": local_index,
                    "gate": float(gates[basis_index]),
                    **params,
                }
            )
        return summaries

    def active_bases(
        self,
        freqs_hz: np.ndarray | list[float] | None = None,
        *,
        threshold: float | None = None,
    ) -> list[YAdmittanceURNActiveBasis]:
        """Return bases selected by output ablation."""

        cutoff = self.config.active_threshold if threshold is None else float(threshold)
        importance = self.output_ablation_importance(freqs_hz)
        summaries = self.parameter_summary()
        active: list[YAdmittanceURNActiveBasis] = []
        for item, imp in zip(summaries, importance, strict=True):
            if imp < cutoff:
                continue
            params = {
                k: float(v)
                for k, v in item.items()
                if k not in {"basis_index", "basis_type", "local_index", "gate"}
            }
            active.append(
                YAdmittanceURNActiveBasis(
                    basis_index=int(item["basis_index"]),
                    basis_type=str(item["basis_type"]),
                    local_index=int(item["local_index"]),
                    importance=float(imp),
                    gate=float(item["gate"]),
                    parameters=params,
                )
            )
        active.sort(key=lambda x: x.importance, reverse=True)
        return active

    def _parameter_values_by_type(self) -> dict[str, list[dict[str, float]]]:
        with torch.no_grad():
            out: dict[str, list[dict[str, float]]] = {
                "debye": [{"tau": float(x)} for x in self._tau(self.debye_tau_raw).cpu().numpy()],
                "magnetic_debye": [
                    {"tau": float(x)} for x in self._tau(self.mag_debye_tau_raw).cpu().numpy()
                ],
                "cole_cole": [],
                "magnetic_cole_cole": [],
                "inductive_cpe": [],
                "capacitive_cpe": [],
                "series_rlc": [],
                "parallel_rlc": [],
                "coil_antiresonance": [],
            }
            cc_tau = self._tau(self.cc_tau_raw).cpu().numpy()
            cc_alpha = _bounded_from_raw(
                self.cc_alpha_raw, self.config.alpha_min, self.config.alpha_max
            ).cpu().numpy()
            out["cole_cole"] = [
                {"tau": float(t), "alpha": float(a)} for t, a in zip(cc_tau, cc_alpha, strict=True)
            ]
            mag_cc_tau = self._tau(self.mag_cc_tau_raw).cpu().numpy()
            mag_cc_alpha = _bounded_from_raw(
                self.mag_cc_alpha_raw, self.config.alpha_min, self.config.alpha_max
            ).cpu().numpy()
            out["magnetic_cole_cole"] = [
                {"tau": float(t), "alpha": float(a)}
                for t, a in zip(mag_cc_tau, mag_cc_alpha, strict=True)
            ]
            out["inductive_cpe"] = [
                {"beta": float(x)}
                for x in _bounded_from_raw(
                    self.ind_cpe_beta_raw, self.config.beta_min, self.config.beta_max
                ).cpu().numpy()
            ]
            out["capacitive_cpe"] = [
                {"beta": float(x)}
                for x in _bounded_from_raw(
                    self.cap_cpe_beta_raw, self.config.beta_min, self.config.beta_max
                ).cpu().numpy()
            ]
            omega0 = self._omega0(self.rlc_omega0_raw).cpu().numpy()
            q = _bounded_from_raw(self.rlc_q_raw, self.config.q_min, self.config.q_max).cpu().numpy()
            out["series_rlc"] = [
                {"omega0": float(w0), "f0": float(w0 / (2.0 * np.pi)), "q": float(q_i)}
                for w0, q_i in zip(omega0, q, strict=True)
            ]
            prlc_omega0 = self._omega0(self.prlc_omega0_raw).cpu().numpy()
            prlc_q = _bounded_from_raw(
                self.prlc_q_raw, self.config.q_min, self.config.q_max
            ).cpu().numpy()
            out["parallel_rlc"] = [
                {"omega0": float(w0), "f0": float(w0 / (2.0 * np.pi)), "q": float(q_i)}
                for w0, q_i in zip(prlc_omega0, prlc_q, strict=True)
            ]
            coil_omega0 = self._omega0(self.coil_omega0_raw).cpu().numpy()
            coil_q = _bounded_from_raw(
                self.coil_q_raw, self.config.q_min, self.config.q_max
            ).cpu().numpy()
            out["coil_antiresonance"] = [
                {"omega0": float(w0), "f0": float(w0 / (2.0 * np.pi)), "q": float(q_i)}
                for w0, q_i in zip(coil_omega0, coil_q, strict=True)
            ]
        return out


def train_y_admittance_urn(
    freqs_hz: np.ndarray | list[float],
    z_data: np.ndarray | list[complex],
    config: YAdmittanceURNConfig | None = None,
    *,
    warm_start_model: YAdmittanceURN | None = None,
    keep_best_epoch: bool = True,
    best_check_interval: int = 1,
    verbose: bool = True,
) -> YAdmittanceURN:
    """Train the SA/RM-2026 Y-domain URN and return the best model."""

    freqs = _as_1d_float64(freqs_hz, "freqs_hz")
    z_arr = _as_1d_complex128(z_data, "z_data")
    if freqs.shape != z_arr.shape:
        raise ValueError("freqs_hz and z_data must have the same length")
    cfg = config or YAdmittanceURNConfig()
    omega = torch.tensor(2.0 * np.pi * freqs, dtype=torch.float64)
    z_target = torch.tensor(z_arr, dtype=torch.complex128)

    best_model: YAdmittanceURN | None = None
    best_loss = float("inf")
    best_history: list[dict[str, float]] = []
    check_interval = max(int(best_check_interval), 1)
    for restart in range(cfg.n_restarts):
        torch.manual_seed(int(cfg.seed) + 7919 * restart)
        model = YAdmittanceURN(freqs, cfg, z_data=z_arr)
        if warm_start_model is not None:
            model.initialize_from_model(warm_start_model)
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
            loss, parts = model.loss(omega, z_target)
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

        final_loss, final_parts = model.loss(omega, z_target)
        final_value = float(final_loss.detach().cpu())
        if verbose:
            active_count = len(model.active_bases(freqs))
            print(
                f"  Y-URN restart {restart + 1}/{cfg.n_restarts}: "
                f"loss={final_value:.6e}, fit={final_parts['loss_fit']:.6e}, "
                f"active={active_count}"
            )
        if final_value < best_loss:
            best_loss = final_value
            best_model = YAdmittanceURN(freqs, cfg, z_data=z_arr)
            best_model.load_state_dict(model.state_dict())
            best_history = history

    if best_model is None:
        raise RuntimeError("Y-admittance URN training did not produce a model")
    best_model.training_history = best_history
    return best_model


def refit_y_admittance_active_bases(
    freqs_hz: np.ndarray | list[float],
    z_data: np.ndarray | list[complex],
    active_bases: list[YAdmittanceURNActiveBasis],
    config: YAdmittanceURNConfig | None = None,
    *,
    verbose: bool = True,
) -> YAdmittanceURN:
    """Refit only selected active bases with frequency-independent gates."""

    freqs = _as_1d_float64(freqs_hz, "freqs_hz")
    z_arr = _as_1d_complex128(z_data, "z_data")
    active = list(active_bases)
    if freqs.shape != z_arr.shape:
        raise ValueError("freqs_hz and z_data must have the same length")

    cfg = active_basis_refit_config(active, config)
    tau_values = [
        float(item.parameters["tau"])
        for item in active
        if "tau" in item.parameters and float(item.parameters["tau"]) > 0.0
    ]
    if tau_values:
        omega_min = 2.0 * np.pi * float(np.min(freqs))
        omega_max = 2.0 * np.pi * float(np.max(freqs))
        auto_tau_min = cfg.tau_min or (0.1 / omega_max)
        auto_tau_max = cfg.tau_max or (10.0 / omega_min)
        cfg.tau_min = min(auto_tau_min, min(tau_values) * 0.25)
        cfg.tau_max = max(auto_tau_max, max(tau_values) * 4.0)

    omega = torch.tensor(2.0 * np.pi * freqs, dtype=torch.float64)
    z_target = torch.tensor(z_arr, dtype=torch.complex128)
    best_model: YAdmittanceURN | None = None
    best_loss = float("inf")
    best_history: list[dict[str, float]] = []
    for restart in range(cfg.n_restarts):
        torch.manual_seed(int(cfg.seed) + 7919 * restart)
        model = YAdmittanceURN(freqs, cfg, z_data=z_arr)
        model.initialize_from_active_bases(active)
        optimizer = optim.Adam(model.parameters(), lr=cfg.lr)
        history: list[dict[str, float]] = []
        for epoch in range(cfg.n_epochs):
            optimizer.zero_grad()
            loss, parts = model.loss(omega, z_target)
            loss.backward()
            torch.nn.utils.clip_grad_norm_(model.parameters(), max_norm=10.0)
            optimizer.step()
            if epoch == 0 or epoch == cfg.n_epochs - 1:
                parts = dict(parts)
                parts["epoch"] = float(epoch + 1)
                history.append(parts)

        final_loss, final_parts = model.loss(omega, z_target)
        final_value = float(final_loss.detach().cpu())
        if verbose:
            print(
                f"  Y-URN active refit {restart + 1}/{cfg.n_restarts}: "
                f"loss={final_value:.6e}, fit={final_parts['loss_fit']:.6e}, "
                f"bases={cfg.total_basis_functions}"
            )
        if final_value < best_loss:
            best_loss = final_value
            best_model = YAdmittanceURN(freqs, cfg, z_data=z_arr)
            best_model.initialize_from_active_bases(active)
            best_model.load_state_dict(model.state_dict())
            best_history = history

    if best_model is None:
        raise RuntimeError("Y-admittance active-basis refit did not produce a model")
    best_model.training_history = best_history
    return best_model


__all__ = [
    "YAdmittanceURN",
    "YAdmittanceURNActiveBasis",
    "YAdmittanceURNConfig",
    "active_basis_refit_config",
    "complex_smooth_l1",
    "refit_y_admittance_active_bases",
    "s_domain_rmse",
    "scattering_transform",
    "train_y_admittance_urn",
]
