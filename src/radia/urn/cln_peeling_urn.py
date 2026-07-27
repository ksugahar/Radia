"""Progressive CLN peeling with paired 22-basis branch identification.

At stage ``n`` the measured driving-point impedance is represented as

``R_n = Z_2n + (Z_2n+1 || R_n+1)``.

A 22-basis open-tail fit supplies the physical basis shapes for the pair.  Each
fitted coefficient is divided between the even series branch and the odd shunt
branch without duplicating the 22 physical basis shapes.  A fresh 22-basis
look-ahead model represents the next tail while the split is selected.  Once
accepted, the two outer branches are frozen and the measured tail is peeled by
the exact inverse map.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Literal

import numpy as np
import torch
from torch import nn, optim

from .y_admittance_urn import (
    _EPS,
    _as_1d_complex128,
    _as_1d_float64,
    complex_smooth_l1,
    scattering_transform,
    s_domain_rmse,
    YAdmittanceURN,
    YAdmittanceURNConfig,
)


Termination = Literal["stored", "lookahead", "constant", "open", "short"]


@dataclass
class CLNPeelingConfig:
    """Configuration for paired-basis CLN peeling.

    The trust-region fields protect the stage-to-stage peeling division.  The
    exact peel ``R_n+1 = [1/(R_n - Z_2n) - 1/Z_2n+1]^-1`` amplifies the
    measurement error of ``R_n`` wherever ``R_n - Z_2n`` cancels (series
    cancellation) or wherever the two admittances nearly cancel so the tail
    barely loads the ladder (parallel-resonance band).  Points inside those
    bands receive a small per-frequency trust weight; the next stage fits and
    is judged on the trusted points only, while the stored exact tail remains
    unmodified for identity reconstruction.
    """

    n_stages: int = 1
    branch_epochs: int = 1200
    branch_lr: float = 4.0e-3
    branch_restarts: int = 1
    branch_sparsity_weight: float = 1.0e-6
    pair_epochs: int = 1600
    pair_polish_epochs: int = 400
    pair_lr: float = 3.0e-3
    pair_restarts: int = 2
    huber_delta: float = 1.0e-2
    gate_init: float = 0.25
    seed: int = 23

    tail_consistency_weight: float = 0.25
    positive_real_weight: float = 1.0
    singularity_weight: float = 2.0e-2
    tail_activity_weight: float = 2.0e-2
    split_binary_weight: float = 0.0
    split_balance_weight: float = 5.0e-2
    seed_regularization_weight: float = 1.0e-4
    min_branch_share: float = 0.08
    min_tail_sensitivity: float = 1.0e-3
    min_tail_admittance_relative: float = 1.0e-4
    residual_real_tolerance: float = 5.0e-2
    hard_split: bool = False
    max_stage_relative_degradation: float = 0.10

    denominator_margin_relative: float = 3.0e-2
    denominator_margin_weight: float = 1.0
    tail_admittance_margin_relative: float = 5.0e-3
    min_trusted_fraction: float = 0.8
    trust_weight_floor: float = 0.0

    # Anti-resonance dictionary extensions (0 = the paper 22-basis dictionary).
    n_parallel_rlc: int = 0
    n_coil_antiresonance: int = 0

    # Warm-start the next stage's composite seed from the previous stage's
    # lookahead branch.  The previous lookahead was optimised through the
    # whole-ladder loss, so starting from it (instead of refitting the peeled
    # tail from scratch) restores the monotonic-improvement property that a
    # greedy frozen-stage ladder otherwise loses.
    composite_warm_start: bool = True

    # Train each stage's pair through the frozen outer branches against the
    # measured response instead of against the peeled local tail.  The frozen
    # branches enter the composition as constants (no gradient flows to past
    # stages), so the freeze principle is untouched, but the optimisation
    # objective is the whole-ladder fit -- a stage can then never improve
    # locally while degrading globally.  The peeled tail remains the seed and
    # the acceptance diagnostics.
    global_objective: bool = True


@dataclass
class CLNBranchFit:
    """Frozen basis-sum impedance for one CLN branch."""

    model: YAdmittanceURN
    response_ref: float
    coefficients: np.ndarray
    role: str
    inverse_sum: bool = True
    final_loss: float = float("nan")

    def response(self, freqs_hz: np.ndarray | list[float]) -> np.ndarray:
        freqs = _as_1d_float64(freqs_hz, "freqs_hz")
        omega = torch.tensor(2.0 * np.pi * freqs, dtype=torch.float64)
        coeff = torch.tensor(np.asarray(self.coefficients), dtype=torch.float64)
        with torch.no_grad():
            basis = self.model.basis_matrix(omega, normalize=True)
            summed = self.response_ref * torch.sum(
                coeff[None, :].to(torch.complex128) * basis,
                dim=1,
            )
            values = 1.0 / (summed + _EPS) if self.inverse_sum else summed
        return values.detach().cpu().numpy()

    def active_count(self, threshold: float = 1.0e-2) -> int:
        coeff = np.asarray(self.coefficients, dtype=np.float64)
        if coeff.size == 0 or np.max(coeff) <= 0.0:
            return 0
        return int(np.count_nonzero(coeff > threshold * np.max(coeff)))


@dataclass
class CLNPeelingStage:
    """One frozen even-series/odd-shunt CLN stage.

    ``sample_weight`` is the per-frequency trust inherited from the previous
    stages (all ones at stage 0); ``tail_trust_weight`` is the trust this
    stage assigns to its own exactly peeled tail, which the next stage
    multiplies into its inherited weight.
    """

    index: int
    composite_seed: CLNBranchFit
    series_impedance: CLNBranchFit
    shunt_impedance: CLNBranchFit
    lookahead_tail: CLNBranchFit
    split_fraction: np.ndarray
    exact_tail_impedance: np.ndarray
    accepted: bool
    sample_weight: np.ndarray | None = None
    tail_trust_weight: np.ndarray | None = None
    metrics: dict[str, float | int | bool] = field(default_factory=dict)


@dataclass
class CLNPeelingNetwork:
    """Progressively identified CLN on its training frequency grid."""

    freqs_hz: np.ndarray
    stages: list[CLNPeelingStage]
    final_tail_impedance: np.ndarray
    training_log: list[dict[str, float | int | bool]] = field(default_factory=list)

    def predict(self, freqs_hz: np.ndarray | list[float]) -> np.ndarray:
        return self.predict_terminated(freqs_hz, termination="stored")

    def predict_terminated(
        self,
        freqs_hz: np.ndarray | list[float],
        *,
        termination: Termination = "stored",
    ) -> np.ndarray:
        """Evaluate the ladder with a stored, learned, or simple tail.

        The frozen branches are basis models, so every termination except
        ``stored`` accepts an arbitrary frequency grid.  ``stored`` uses the
        exactly peeled measurement table and is therefore only defined on the
        training grid; it exists for identity-reconstruction checks and must
        not be reported as a fit result.
        """

        freqs = _as_1d_float64(freqs_hz, "freqs_hz")
        on_training_grid = freqs.shape == self.freqs_hz.shape and np.allclose(
            freqs, self.freqs_hz
        )
        if termination == "stored":
            if not on_training_grid:
                raise ValueError(
                    "termination='stored' evaluates the exactly peeled tail table "
                    "and is only defined on the training grid"
                )
            z_tail = np.asarray(self.final_tail_impedance, dtype=np.complex128)
        elif termination == "lookahead":
            if self.stages:
                z_tail = self.stages[-1].lookahead_tail.response(freqs)
            else:
                if not on_training_grid:
                    raise ValueError(
                        "termination='lookahead' without accepted stages falls back "
                        "to the stored tail table and needs the training grid"
                    )
                z_tail = np.asarray(self.final_tail_impedance, dtype=np.complex128)
        elif termination == "constant":
            z_tail = np.full(
                freqs.shape,
                np.median(np.asarray(self.final_tail_impedance, dtype=np.complex128)),
                dtype=np.complex128,
            )
        elif termination == "open":
            z_tail = np.full(freqs.shape, 1.0e30 + 0.0j, dtype=np.complex128)
        elif termination == "short":
            z_tail = np.zeros(freqs.shape, dtype=np.complex128)
        else:
            raise ValueError("termination must be stored, lookahead, constant, open, or short")

        for stage in reversed(self.stages):
            z_tail = _compose_stage(
                stage.series_impedance.response(freqs),
                stage.shunt_impedance.response(freqs),
                z_tail,
            )
        return z_tail

    def audit_passivity(
        self,
        freqs_hz: np.ndarray | list[float] | None = None,
        *,
        points_per_decade: int = 64,
        extrapolation_decades: float = 1.0,
        tolerance: float = 0.0,
    ) -> dict:
        """Positive-real audit of every frozen branch on a dense grid.

        Evaluates each frozen series/shunt/lookahead branch and the full
        lookahead ladder on a dense logarithmic grid (default: denser than the
        training grid and extended by ``extrapolation_decades`` on both sides)
        and reports minimum real parts, their frequencies, and finiteness.
        The stored tail table is excluded on purpose: it is the exactly peeled
        measurement residue, not a physical circuit.
        """

        if freqs_hz is None:
            lo = float(np.min(self.freqs_hz)) / (10.0 ** float(extrapolation_decades))
            hi = float(np.max(self.freqs_hz)) * (10.0 ** float(extrapolation_decades))
            n = max(int(np.ceil(np.log10(hi / lo) * points_per_decade)) + 1, 2)
            freqs = np.logspace(np.log10(lo), np.log10(hi), n)
        else:
            freqs = _as_1d_float64(freqs_hz, "freqs_hz")

        def branch_report(stage_index: int, branch: CLNBranchFit) -> dict:
            values = branch.response(freqs)
            admittance = 1.0 / (values + _EPS)
            finite = bool(
                np.all(np.isfinite(values.real)) and np.all(np.isfinite(values.imag))
            )
            return {
                "stage": stage_index,
                "role": branch.role,
                "finite": finite,
                "min_re_z": float(np.min(values.real)),
                "min_re_z_freq_hz": float(freqs[int(np.argmin(values.real))]),
                "min_re_y": float(np.min(admittance.real)),
                "min_re_y_freq_hz": float(freqs[int(np.argmin(admittance.real))]),
            }

        branches = []
        for stage in self.stages:
            branches.append(branch_report(stage.index, stage.series_impedance))
            branches.append(branch_report(stage.index, stage.shunt_impedance))
            branches.append(branch_report(stage.index, stage.lookahead_tail))

        report: dict = {
            "freq_min_hz": float(freqs[0]),
            "freq_max_hz": float(freqs[-1]),
            "n_points": int(freqs.size),
            "tolerance": float(tolerance),
            "branches": branches,
        }
        checks = [
            b["finite"] and b["min_re_z"] >= -tolerance and b["min_re_y"] >= -tolerance
            for b in branches
        ]
        if self.stages:
            ladder = self.predict_terminated(freqs, termination="lookahead")
            ladder_report = {
                "finite": bool(
                    np.all(np.isfinite(ladder.real)) and np.all(np.isfinite(ladder.imag))
                ),
                "min_re_z": float(np.min(ladder.real)),
                "min_re_z_freq_hz": float(freqs[int(np.argmin(ladder.real))]),
            }
            report["ladder_lookahead"] = ladder_report
            checks.append(
                ladder_report["finite"] and ladder_report["min_re_z"] >= -tolerance
            )
        report["passive"] = bool(all(checks)) if checks else True
        return report

    def s_domain_rmse(self, z_target: np.ndarray | list[complex]) -> float:
        return s_domain_rmse(self.predict(self.freqs_hz), z_target)

    def s_domain_rmse_terminated(
        self,
        z_target: np.ndarray | list[complex],
        *,
        termination: Termination = "lookahead",
    ) -> float:
        return s_domain_rmse(
            self.predict_terminated(self.freqs_hz, termination=termination),
            z_target,
        )


def _branch_config(config: CLNPeelingConfig) -> YAdmittanceURNConfig:
    return YAdmittanceURNConfig.paper_22_basis(
        n_epochs=config.branch_epochs,
        lr=config.branch_lr,
        n_restarts=1,
        sparsity_weight=0.0,
        gate_init=config.gate_init,
        huber_delta=config.huber_delta,
        n_parallel_rlc=config.n_parallel_rlc,
        n_coil_antiresonance=config.n_coil_antiresonance,
    )


def _direct_basis_response(
    model: YAdmittanceURN,
    omega: torch.Tensor,
    response_ref: float,
    coefficients: torch.Tensor,
    *,
    inverse_sum: bool = True,
) -> torch.Tensor:
    basis = model.basis_matrix(omega, normalize=True)
    summed = response_ref * torch.sum(
        coefficients[None, :].to(torch.complex128) * basis,
        dim=1,
    )
    return 1.0 / (summed + _EPS) if inverse_sum else summed


def _s_loss(
    pred: torch.Tensor,
    target: torch.Tensor,
    z0: float,
    huber_delta: float,
    weight: torch.Tensor | None = None,
) -> torch.Tensor:
    residual = scattering_transform(pred, z0) - scattering_transform(target, z0)
    return complex_smooth_l1(residual, huber_delta, weight=weight)


def _weighted_mean(values: torch.Tensor, weight: torch.Tensor | None) -> torch.Tensor:
    if weight is None:
        return torch.mean(values)
    w = weight.to(dtype=values.dtype, device=values.device)
    return torch.sum(w * values) / torch.sum(w).clamp_min(_EPS)


def _np_weighted_mean(values: np.ndarray, weight: np.ndarray | None) -> float:
    if weight is None:
        return float(np.mean(values))
    w = np.asarray(weight, dtype=np.float64)
    return float(np.sum(w * values) / max(float(np.sum(w)), _EPS))


def _weighted_s_domain_rmse(
    z_fit: np.ndarray,
    z_target: np.ndarray,
    weight: np.ndarray | None,
    *,
    z0: float | None = None,
) -> float:
    """S-domain RMSE with nonnegative per-frequency weights."""

    if weight is None:
        return s_domain_rmse(z_fit, z_target, z0=z0)
    fit = np.asarray(z_fit, dtype=np.complex128)
    target = np.asarray(z_target, dtype=np.complex128)
    z0_value = float(np.median(np.abs(target))) if z0 is None else float(z0)
    z0_value = max(z0_value, _EPS)
    s_fit = (fit - z0_value) / (fit + z0_value + _EPS)
    s_target = (target - z0_value) / (target + z0_value + _EPS)
    return float(
        np.sqrt(_np_weighted_mean(np.abs(s_fit - s_target) ** 2, weight))
    )


def _tail_trust_weight(
    target: np.ndarray,
    series: np.ndarray,
    tail_admittance: np.ndarray,
    z_ref: float,
    config: CLNPeelingConfig,
) -> np.ndarray:
    """Per-frequency confidence of the exactly peeled tail.

    Two error-amplifying mechanisms make individual peeled-tail values
    meaningless:

    1. ``R_n - Z_2n`` cancels, so ``1/(R_n - Z_2n)`` amplifies the measurement
       error of ``R_n`` (series cancellation).
    2. ``1/(R_n - Z_2n)`` and ``1/Z_2n+1`` nearly cancel, so the peeled tail
       admittance is tiny: the tail barely loads the ladder there and its
       peeled value is amplified noise with an unstable sign (the
       parallel-resonance band that produced the PCB -59.7 kOhm spike).

    Both are mapped to smooth weights in ``(0, 1]`` and multiplied.
    """

    parallel_part = target - series
    cancellation = np.abs(parallel_part) / (
        np.abs(target) + np.abs(series) + _EPS
    )
    m_series = float(config.denominator_margin_relative)
    w_series = cancellation**2 / (cancellation**2 + m_series**2)

    activity = np.abs(tail_admittance) * float(z_ref)
    m_tail = float(config.tail_admittance_margin_relative)
    w_tail = activity**2 / (activity**2 + m_tail**2)

    floor = float(config.trust_weight_floor)
    return floor + (1.0 - floor) * (w_series * w_tail)


def _inverse_softplus(values: torch.Tensor) -> torch.Tensor:
    values = values.clamp_min(1.0e-12)
    return torch.log(torch.expm1(values).clamp_min(1.0e-30))


def _compose_stage(
    series: np.ndarray | torch.Tensor,
    shunt: np.ndarray | torch.Tensor,
    tail: np.ndarray | torch.Tensor,
) -> np.ndarray | torch.Tensor:
    """Compose ``series + (shunt || tail)`` without changing array type."""

    return series + 1.0 / (1.0 / (shunt + _EPS) + 1.0 / (tail + _EPS) + _EPS)


def _peel_tail(
    driving_point: np.ndarray | torch.Tensor,
    series: np.ndarray | torch.Tensor,
    shunt: np.ndarray | torch.Tensor,
) -> tuple[np.ndarray | torch.Tensor, np.ndarray | torch.Tensor]:
    """Return exact tail impedance and tail admittance for a frozen pair."""

    parallel_part = driving_point - series
    tail_admittance = 1.0 / (parallel_part + _EPS) - 1.0 / (shunt + _EPS)
    return 1.0 / (tail_admittance + _EPS), tail_admittance


def _fit_composite_seed(
    freqs_hz: np.ndarray,
    target: np.ndarray,
    config: CLNPeelingConfig,
    *,
    seed_offset: int,
    sample_weight: np.ndarray | None = None,
    warm_start_model: YAdmittanceURN | None = None,
) -> CLNBranchFit:
    target_arr = _as_1d_complex128(target, "target")
    omega = torch.tensor(2.0 * np.pi * freqs_hz, dtype=torch.float64)
    target_t = torch.tensor(target_arr, dtype=torch.complex128)
    weight_t = (
        None
        if sample_weight is None
        else torch.tensor(np.asarray(sample_weight, dtype=np.float64))
    )
    best_loss = float("inf")
    best_model: YAdmittanceURN | None = None
    cfg = _branch_config(config)

    for restart in range(config.branch_restarts):
        torch.manual_seed(int(config.seed) + seed_offset + 7919 * restart)
        model = YAdmittanceURN(freqs_hz, cfg, z_data=target_arr)
        if warm_start_model is not None and restart == 0:
            # Reproduce the previous lookahead branch exactly: copy its
            # parameters and rescale the gates for the new y_ref, so
            # y_ref_new * sum(g_new B) == y_ref_old * sum(g_old B).
            model.load_state_dict(warm_start_model.state_dict())
            scale = float(warm_start_model.y_ref) / float(model.y_ref)
            with torch.no_grad():
                gates = torch.nn.functional.softplus(model.gate_raw) * scale
                model.gate_raw.copy_(
                    torch.log(torch.expm1(gates.clamp_min(1.0e-12)))
                )
        optimizer = optim.Adam(model.parameters(), lr=config.branch_lr)
        for _epoch in range(config.branch_epochs):
            optimizer.zero_grad()
            pred = _direct_basis_response(
                model,
                omega,
                model.y_ref,
                model.gates(),
                inverse_sum=True,
            )
            loss_fit = _s_loss(
                pred, target_t, model.z0, config.huber_delta, weight=weight_t
            )
            loss_sparse = config.branch_sparsity_weight * torch.mean(model.gates())
            loss = loss_fit + loss_sparse
            loss.backward()
            torch.nn.utils.clip_grad_norm_(model.parameters(), max_norm=10.0)
            optimizer.step()
        with torch.no_grad():
            pred = _direct_basis_response(
                model,
                omega,
                model.y_ref,
                model.gates(),
                inverse_sum=True,
            )
            loss = float(
                _s_loss(pred, target_t, model.z0, config.huber_delta, weight=weight_t)
                .detach()
                .cpu()
            )
        if loss < best_loss:
            best_loss = loss
            best_model = YAdmittanceURN(freqs_hz, cfg, z_data=target_arr)
            best_model.load_state_dict(model.state_dict())

    if best_model is None:
        raise RuntimeError("CLN composite fit did not produce a model")
    return CLNBranchFit(
        model=best_model,
        response_ref=best_model.y_ref,
        coefficients=best_model.gates().detach().cpu().numpy(),
        role="composite_seed",
        inverse_sum=True,
        final_loss=best_loss,
    )


def _split_fraction(
    split_logits: torch.Tensor,
    *,
    hard: bool,
) -> torch.Tensor:
    soft = torch.sigmoid(split_logits)
    if not hard:
        return soft
    hard_value = (soft >= 0.5).to(torch.float64)
    return hard_value + soft - soft.detach()


def _ensure_nonempty_mask(mask: np.ndarray, logits: np.ndarray) -> np.ndarray:
    result = np.asarray(mask, dtype=bool).copy()
    if not np.any(result):
        result[int(np.argmax(logits))] = True
    if np.all(result):
        result[int(np.argmin(logits))] = False
    return result


def _pair_loss(
    *,
    target: torch.Tensor,
    outer_basis: torch.Tensor,
    total_coefficients: torch.Tensor,
    split_fraction: torch.Tensor,
    child_response: torch.Tensor,
    outer_admittance_ref: float,
    impedance_ref: float,
    seed_coefficients: torch.Tensor,
    config: CLNPeelingConfig,
    sample_weight: torch.Tensor | None = None,
    outer_chain: list[tuple[torch.Tensor, torch.Tensor]] | None = None,
    measured: torch.Tensor | None = None,
    measured_ref: float | None = None,
) -> tuple[torch.Tensor, dict[str, torch.Tensor]]:
    even_coeff = total_coefficients * split_fraction
    odd_coeff = total_coefficients * (1.0 - split_fraction)
    series_y = outer_admittance_ref * torch.sum(
        even_coeff[None, :].to(torch.complex128) * outer_basis,
        dim=1,
    )
    shunt_y = outer_admittance_ref * torch.sum(
        odd_coeff[None, :].to(torch.complex128) * outer_basis,
        dim=1,
    )
    series = 1.0 / (series_y + _EPS)
    shunt = 1.0 / (shunt_y + _EPS)
    pred = _compose_stage(series, shunt, child_response)
    exact_tail, exact_tail_y = _peel_tail(target, series, shunt)
    parallel_part = target - series

    if outer_chain is not None and measured is not None and measured_ref is not None:
        # Whole-ladder objective: compose the current pair through the frozen
        # outer branches (constants -- no gradient reaches past stages) and
        # fit the raw measured response on every point.
        pred_global = pred
        for outer_series, outer_shunt in reversed(outer_chain):
            pred_global = _compose_stage(outer_series, outer_shunt, pred_global)
        loss_fit = _s_loss(pred_global, measured, measured_ref, config.huber_delta)
    else:
        loss_fit = _s_loss(
            pred, target, impedance_ref, config.huber_delta, weight=sample_weight
        )
    loss_tail = _s_loss(
        child_response, exact_tail, impedance_ref, config.huber_delta, weight=sample_weight
    )

    tol = float(config.residual_real_tolerance)
    parallel_real_norm = parallel_part.real / impedance_ref
    tail_y_real_norm = exact_tail_y.real * impedance_ref
    # Weighted mean penalties only.  A worst-point (max) penalty was tried on
    # 2026-07-27 to close the gap between the continuous training weights and
    # the binary trusted-min acceptance gate, but on the PCB stage-2 data the
    # concentrated gradient destroyed the fit (parent RMSE 8.5e-3 -> 7.7e-2)
    # while the trusted violation worsened (-1.14 -> -6.4); do not reintroduce.
    pr_parallel = _weighted_mean(
        torch.relu(-parallel_real_norm - tol) ** 2, sample_weight
    )
    pr_tail_y = _weighted_mean(torch.relu(-tail_y_real_norm - tol) ** 2, sample_weight)

    min_abs = 1.0e-5
    singular_parallel = torch.mean(
        torch.relu(min_abs - torch.abs(parallel_part) / impedance_ref) ** 2
    )
    singular_shunt = torch.mean(
        torch.relu(min_abs - torch.abs(shunt) / impedance_ref) ** 2
    )
    cancellation = torch.abs(parallel_part) / (
        torch.abs(target) + torch.abs(series) + _EPS
    )
    denominator_margin = _weighted_mean(
        torch.relu(
            torch.as_tensor(
                config.denominator_margin_relative, dtype=torch.float64
            )
            - cancellation
        )
        ** 2,
        sample_weight,
    )
    tail_activity = torch.mean(torch.abs(exact_tail_y) * impedance_ref)
    inactive_tail = torch.relu(
        torch.as_tensor(config.min_tail_admittance_relative, dtype=torch.float64)
        - tail_activity
    ) ** 2

    coefficient_sum = torch.sum(total_coefficients).clamp_min(1.0e-24)
    series_share = torch.sum(total_coefficients * split_fraction) / coefficient_sum
    min_share = torch.as_tensor(config.min_branch_share, dtype=torch.float64)
    balance = torch.relu(min_share - series_share) ** 2 + torch.relu(
        series_share - (1.0 - min_share)
    ) ** 2
    binary = torch.mean(split_fraction * (1.0 - split_fraction))
    seed_reg = torch.mean(
        (
            torch.log(total_coefficients.clamp_min(1.0e-12))
            - torch.log(seed_coefficients.clamp_min(1.0e-12))
        )
        ** 2
    )

    loss = (
        loss_fit
        + config.tail_consistency_weight * loss_tail
        + config.positive_real_weight * (pr_parallel + pr_tail_y)
        + config.singularity_weight * (singular_parallel + singular_shunt)
        + config.denominator_margin_weight * denominator_margin
        + config.tail_activity_weight * inactive_tail
        + config.split_binary_weight * binary
        + config.split_balance_weight * balance
        + config.seed_regularization_weight * seed_reg
    )
    return loss, {
        "series": series,
        "shunt": shunt,
        "pred": pred,
        "exact_tail": exact_tail,
        "exact_tail_y": exact_tail_y,
        "loss_fit": loss_fit,
        "loss_tail": loss_tail,
        "series_share": series_share,
        "tail_activity": tail_activity,
    }


def _fit_paired_stage(
    freqs_hz: np.ndarray,
    target: np.ndarray,
    composite_seed: CLNBranchFit,
    config: CLNPeelingConfig,
    *,
    stage_index: int,
    sample_weight: np.ndarray | None = None,
    outer_chain: list[tuple[np.ndarray, np.ndarray]] | None = None,
    measured: np.ndarray | None = None,
    measured_ref: float | None = None,
) -> CLNPeelingStage:
    omega = torch.tensor(2.0 * np.pi * freqs_hz, dtype=torch.float64)
    target_t = torch.tensor(target, dtype=torch.complex128)
    outer_admittance_ref = composite_seed.response_ref
    impedance_ref = composite_seed.model.z0
    seed_coeff = torch.tensor(composite_seed.coefficients, dtype=torch.float64)
    sample_weight_np = (
        None
        if sample_weight is None
        else np.asarray(sample_weight, dtype=np.float64)
    )
    weight_t = (
        None if sample_weight_np is None else torch.tensor(sample_weight_np)
    )
    outer_chain_t = (
        None
        if outer_chain is None
        else [
            (
                torch.tensor(np.asarray(series, dtype=np.complex128)),
                torch.tensor(np.asarray(shunt, dtype=np.complex128)),
            )
            for series, shunt in outer_chain
        ]
    )
    measured_t = (
        None
        if measured is None
        else torch.tensor(np.asarray(measured, dtype=np.complex128))
    )
    with torch.no_grad():
        outer_basis = composite_seed.model.basis_matrix(omega, normalize=True).detach()

    best: tuple[float, CLNPeelingStage] | None = None
    cfg = _branch_config(config)
    for restart in range(config.pair_restarts):
        torch.manual_seed(int(config.seed) + 100_000 * stage_index + 7919 * restart)
        child_model = YAdmittanceURN(freqs_hz, cfg, z_data=target)
        child_model.load_state_dict(composite_seed.model.state_dict())
        with torch.no_grad():
            child_model.gate_raw.copy_(_inverse_softplus(0.75 * seed_coeff))

        # This initialization exactly reproduces the open-tail composite B:
        # Z_even=0.2B, Z_odd=2B, and R_next=4B/3.
        total_raw = nn.Parameter(_inverse_softplus(5.5 * seed_coeff).detach().clone())
        initial_fraction = torch.full_like(seed_coeff, 5.0 / 5.5)
        initial_logits = torch.log(initial_fraction / (1.0 - initial_fraction))
        split_logits = nn.Parameter(initial_logits + 0.02 * torch.randn_like(seed_coeff))
        optimizer = optim.Adam(
            [total_raw, split_logits, *child_model.parameters()],
            lr=config.pair_lr,
        )

        for epoch in range(config.pair_epochs):
            optimizer.zero_grad()
            total = torch.nn.functional.softplus(total_raw)
            fraction = _split_fraction(split_logits, hard=config.hard_split)
            child = _direct_basis_response(
                child_model,
                omega,
                child_model.y_ref,
                child_model.gates(),
                inverse_sum=True,
            )
            loss, _parts = _pair_loss(
                target=target_t,
                outer_basis=outer_basis,
                total_coefficients=total,
                split_fraction=fraction,
                child_response=child,
                outer_admittance_ref=outer_admittance_ref,
                impedance_ref=impedance_ref,
                seed_coefficients=seed_coeff,
                config=config,
                sample_weight=weight_t,
                outer_chain=outer_chain_t,
                measured=measured_t,
                measured_ref=measured_ref,
            )
            loss.backward()
            torch.nn.utils.clip_grad_norm_(
                [total_raw, split_logits, *child_model.parameters()],
                max_norm=10.0,
            )
            optimizer.step()

        logits_np = split_logits.detach().cpu().numpy()
        if config.hard_split:
            mask_np = _ensure_nonempty_mask(logits_np >= 0.0, logits_np)
            frozen_fraction_np = mask_np.astype(np.float64)
        else:
            frozen_fraction_np = 1.0 / (1.0 + np.exp(-logits_np))
        frozen_fraction_t = torch.tensor(frozen_fraction_np, dtype=torch.float64)

        polish_optimizer = optim.Adam(
            [total_raw, *child_model.parameters()],
            lr=0.5 * config.pair_lr,
        )
        for _epoch in range(config.pair_polish_epochs):
            polish_optimizer.zero_grad()
            total = torch.nn.functional.softplus(total_raw)
            child = _direct_basis_response(
                child_model,
                omega,
                child_model.y_ref,
                child_model.gates(),
                inverse_sum=True,
            )
            loss, _parts = _pair_loss(
                target=target_t,
                outer_basis=outer_basis,
                total_coefficients=total,
                split_fraction=frozen_fraction_t,
                child_response=child,
                outer_admittance_ref=outer_admittance_ref,
                impedance_ref=impedance_ref,
                seed_coefficients=seed_coeff,
                config=config,
                sample_weight=weight_t,
                outer_chain=outer_chain_t,
                measured=measured_t,
                measured_ref=measured_ref,
            )
            loss.backward()
            torch.nn.utils.clip_grad_norm_(
                [total_raw, *child_model.parameters()],
                max_norm=10.0,
            )
            polish_optimizer.step()

        with torch.no_grad():
            total = torch.nn.functional.softplus(total_raw)
            child = _direct_basis_response(
                child_model,
                omega,
                child_model.y_ref,
                child_model.gates(),
                inverse_sum=True,
            )
            final_loss, parts = _pair_loss(
                target=target_t,
                outer_basis=outer_basis,
                total_coefficients=total,
                split_fraction=frozen_fraction_t,
                child_response=child,
                outer_admittance_ref=outer_admittance_ref,
                impedance_ref=impedance_ref,
                seed_coefficients=seed_coeff,
                config=config,
                sample_weight=weight_t,
                outer_chain=outer_chain_t,
                measured=measured_t,
                measured_ref=measured_ref,
            )

        total_np = total.detach().cpu().numpy()
        even_coeff = total_np * frozen_fraction_np
        odd_coeff = total_np * (1.0 - frozen_fraction_np)
        series_np = parts["series"].detach().cpu().numpy()
        shunt_np = parts["shunt"].detach().cpu().numpy()
        child_np = parts["pred"].detach().cpu().numpy()
        exact_tail_np = parts["exact_tail"].detach().cpu().numpy()
        exact_tail_y_np = parts["exact_tail_y"].detach().cpu().numpy()
        lookahead_np = child.detach().cpu().numpy()

        seed_pred = composite_seed.response(freqs_hz)
        seed_rmse = s_domain_rmse(seed_pred, target)
        parent_rmse = s_domain_rmse(child_np, target)
        seed_rmse_trusted = _weighted_s_domain_rmse(seed_pred, target, sample_weight_np)
        parent_rmse_trusted = _weighted_s_domain_rmse(child_np, target, sample_weight_np)
        tail_rmse = _weighted_s_domain_rmse(
            lookahead_np, exact_tail_np, sample_weight_np, z0=impedance_ref
        )
        parallel_part = target - series_np
        trusted_mask = (
            np.ones(freqs_hz.shape, dtype=bool)
            if sample_weight_np is None
            else sample_weight_np >= 0.5
        )
        if not np.any(trusted_mask):
            trusted_mask = np.ones(freqs_hz.shape, dtype=bool)
        parallel_real_norm = parallel_part.real / impedance_ref
        tail_y_real_norm = exact_tail_y_np.real * impedance_ref
        min_parallel_real = float(np.min(parallel_real_norm))
        min_tail_y_real = float(np.min(tail_y_real_norm))
        trusted_indices = np.flatnonzero(trusted_mask)
        parallel_argmin = trusted_indices[
            int(np.argmin(parallel_real_norm[trusted_mask]))
        ]
        min_parallel_real_trusted = float(parallel_real_norm[parallel_argmin])
        min_parallel_real_trusted_freq = float(freqs_hz[parallel_argmin])
        min_tail_y_real_trusted = float(np.min(tail_y_real_norm[trusted_mask]))
        sensitivity = _np_weighted_mean(
            np.abs(shunt_np / (shunt_np + lookahead_np + _EPS)) ** 2,
            sample_weight_np,
        )
        trust = _tail_trust_weight(
            target, series_np, exact_tail_y_np, impedance_ref, config
        )
        trusted_fraction = _np_weighted_mean(trust, sample_weight_np)
        finite = bool(
            np.all(np.isfinite(exact_tail_np.real))
            and np.all(np.isfinite(exact_tail_np.imag))
        )
        accepted = bool(
            finite
            and min_parallel_real_trusted >= -config.residual_real_tolerance
            and min_tail_y_real_trusted >= -config.residual_real_tolerance
            and sensitivity >= config.min_tail_sensitivity
            and parent_rmse_trusted
            <= seed_rmse_trusted * (1.0 + config.max_stage_relative_degradation)
            + 1.0e-12
            and trusted_fraction >= config.min_trusted_fraction
        )

        child_copy = YAdmittanceURN(freqs_hz, cfg, z_data=target)
        child_copy.load_state_dict(child_model.state_dict())
        series_branch = CLNBranchFit(
            model=composite_seed.model,
            response_ref=outer_admittance_ref,
            coefficients=even_coeff,
            role=f"Z_{2 * stage_index}",
            inverse_sum=True,
            final_loss=float(final_loss.detach().cpu()),
        )
        shunt_branch = CLNBranchFit(
            model=composite_seed.model,
            response_ref=outer_admittance_ref,
            coefficients=odd_coeff,
            role=f"Z_{2 * stage_index + 1}",
            inverse_sum=True,
            final_loss=float(final_loss.detach().cpu()),
        )
        lookahead_branch = CLNBranchFit(
            model=child_copy,
            response_ref=child_copy.y_ref,
            coefficients=child_copy.gates().detach().cpu().numpy(),
            role=f"R_{stage_index + 1}_lookahead",
            inverse_sum=True,
            final_loss=float(parts["loss_tail"].detach().cpu()),
        )
        metrics: dict[str, float | int | bool] = {
            "seed_s_rmse": seed_rmse,
            "parent_lookahead_s_rmse": parent_rmse,
            "seed_s_rmse_trusted": seed_rmse_trusted,
            "parent_lookahead_s_rmse_trusted": parent_rmse_trusted,
            "tail_consistency_s_rmse": tail_rmse,
            "min_parallel_real_normalized": min_parallel_real,
            "min_tail_admittance_real_normalized": min_tail_y_real,
            "min_parallel_real_trusted": min_parallel_real_trusted,
            "min_parallel_real_trusted_freq_hz": min_parallel_real_trusted_freq,
            "min_tail_admittance_real_trusted": min_tail_y_real_trusted,
            "trusted_fraction": trusted_fraction,
            "tail_sensitivity": sensitivity,
            "series_share": float(parts["series_share"].detach().cpu()),
            "tail_activity": float(parts["tail_activity"].detach().cpu()),
            "series_active_count": series_branch.active_count(),
            "shunt_active_count": shunt_branch.active_count(),
            "lookahead_active_count": lookahead_branch.active_count(),
            "accepted": accepted,
        }
        stage = CLNPeelingStage(
            index=stage_index,
            composite_seed=composite_seed,
            series_impedance=series_branch,
            shunt_impedance=shunt_branch,
            lookahead_tail=lookahead_branch,
            split_fraction=frozen_fraction_np,
            exact_tail_impedance=exact_tail_np,
            accepted=accepted,
            sample_weight=(
                np.ones(freqs_hz.shape, dtype=np.float64)
                if sample_weight_np is None
                else sample_weight_np.copy()
            ),
            tail_trust_weight=trust,
            metrics=metrics,
        )
        violation = max(0.0, -min_parallel_real_trusted) + max(
            0.0, -min_tail_y_real_trusted
        )
        score = (
            parent_rmse_trusted
            + config.tail_consistency_weight * tail_rmse
            + 10.0 * violation
        )
        if not accepted:
            score += 1.0
        if best is None or score < best[0]:
            best = (score, stage)

    if best is None:
        raise RuntimeError("CLN paired split did not produce a candidate")
    return best[1]


def train_cln_peeling_urn(
    freqs_hz: np.ndarray | list[float],
    z_data: np.ndarray | list[complex],
    config: CLNPeelingConfig | None = None,
    *,
    verbose: bool = True,
) -> CLNPeelingNetwork:
    """Fit and freeze paired CLN stages from the outside toward the tail."""

    cfg = config or CLNPeelingConfig()
    freqs = _as_1d_float64(freqs_hz, "freqs_hz")
    z_current = _as_1d_complex128(z_data, "z_data")
    if freqs.shape != z_current.shape:
        raise ValueError("freqs_hz and z_data must have the same length")
    if cfg.n_stages <= 0:
        raise ValueError("n_stages must be positive")

    stages: list[CLNPeelingStage] = []
    log: list[dict[str, float | int | bool]] = []
    sample_weight = np.ones(freqs.shape, dtype=np.float64)
    z_measured = z_current.copy()
    measured_ref = float(max(np.median(np.abs(z_measured)), _EPS))
    outer_chain: list[tuple[np.ndarray, np.ndarray]] = []
    previous_global_rmse: float | None = None
    for stage_index in range(cfg.n_stages):
        warm_start = (
            stages[-1].lookahead_tail.model
            if cfg.composite_warm_start and stages
            else None
        )
        composite = _fit_composite_seed(
            freqs,
            z_current,
            cfg,
            seed_offset=10_000 * stage_index,
            sample_weight=sample_weight,
            warm_start_model=warm_start,
        )
        stage = _fit_paired_stage(
            freqs,
            z_current,
            composite,
            cfg,
            stage_index=stage_index,
            sample_weight=sample_weight,
            outer_chain=outer_chain if cfg.global_objective else None,
            measured=z_measured if cfg.global_objective else None,
            measured_ref=measured_ref if cfg.global_objective else None,
        )
        entry: dict[str, float | int | bool] = {
            "stage": stage_index,
            **stage.metrics,
        }
        log.append(entry)
        if verbose:
            print(
                f"  CLN pair {stage_index + 1}/{cfg.n_stages}: "
                f"accepted={stage.accepted}, "
                f"seed={stage.metrics['seed_s_rmse']:.3e}, "
                f"lookahead={stage.metrics['parent_lookahead_s_rmse']:.3e}, "
                f"tail={stage.metrics['tail_consistency_s_rmse']:.3e}, "
                f"trusted={stage.metrics['trusted_fraction']:.3f}, "
                f"Zeven={stage.metrics['series_active_count']}, "
                f"Zodd={stage.metrics['shunt_active_count']}"
            )
        if not stage.accepted:
            break
        # Global stopping rule: a stage replaces the previous lookahead branch
        # by two frozen branches plus a fresh lookahead.  If that replacement
        # degrades the whole-ladder lookahead fit of the measured response,
        # the deeper stage adds no value -- do not freeze it.
        tentative = CLNPeelingNetwork(
            freqs_hz=freqs,
            stages=[*stages, stage],
            final_tail_impedance=np.asarray(
                stage.exact_tail_impedance, dtype=np.complex128
            ),
        )
        global_rmse = s_domain_rmse(
            tentative.predict_terminated(freqs, termination="lookahead"),
            z_measured,
        )
        stage.metrics["global_lookahead_s_rmse"] = global_rmse
        entry["global_lookahead_s_rmse"] = global_rmse
        if (
            previous_global_rmse is not None
            and global_rmse > previous_global_rmse + 1.0e-12
        ):
            stage.metrics["rejected_by_global_degradation"] = True
            entry["rejected_by_global_degradation"] = True
            if verbose:
                print(
                    f"  CLN pair {stage_index + 1}: rejected by global "
                    f"degradation ({global_rmse:.3e} > {previous_global_rmse:.3e})"
                )
            break
        previous_global_rmse = global_rmse
        stages.append(stage)
        outer_chain.append(
            (
                stage.series_impedance.response(freqs),
                stage.shunt_impedance.response(freqs),
            )
        )
        z_current = np.asarray(stage.exact_tail_impedance, dtype=np.complex128)
        sample_weight = sample_weight * np.asarray(
            stage.tail_trust_weight, dtype=np.float64
        )

    return CLNPeelingNetwork(
        freqs_hz=freqs,
        stages=stages,
        final_tail_impedance=z_current,
        training_log=log,
    )


__all__ = [
    "CLNBranchFit",
    "CLNPeelingConfig",
    "CLNPeelingNetwork",
    "CLNPeelingStage",
    "train_cln_peeling_urn",
]
