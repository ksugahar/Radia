"""Measurement-budgeted, duplicate-aware reduction of a fitted single Y-URN.

The budget is caller-owned uncertainty in ohms, not a universal fit target.
This is a greedy candidate path, not a globally minimal circuit certificate.
"""

from dataclasses import replace

import numpy as np
import torch

from .y_admittance_urn import (
    YAdmittanceURN,
    log_component_rmse,
    refit_y_admittance_active_bases,
    s_domain_rmse,
)


def redundant_y_bases(model, freqs_hz, *, response_tolerance=0.01):
    """Find same-family shapes equal up to positive scale on the supplied band.

    Near time constants alone are insufficient: exponent and resonance Q affect
    the full complex response too. Pairs are diagnostics, not physical identity.
    """
    if not np.isfinite(response_tolerance) or not 0 <= response_tolerance < 1:
        raise ValueError("response_tolerance must be finite and in [0, 1)")
    freqs = np.asarray(freqs_hz, dtype=float)
    if (
        freqs.ndim != 1
        or not len(freqs)
        or not np.all(np.isfinite(freqs) & (freqs > 0))
    ):
        raise ValueError("freqs_hz must be a finite positive vector")
    with torch.no_grad():
        columns = (
            model.basis_matrix(torch.tensor(2 * np.pi * freqs, dtype=torch.float64))
            .cpu()
            .numpy()
        )
    norms = np.linalg.norm(columns, axis=0)
    pairs = []
    for i, (family, _) in enumerate(model.basis_labels):
        for j in range(i + 1, len(model.basis_labels)):
            if family != model.basis_labels[j][0] or min(norms[i], norms[j]) <= 0:
                continue
            distance = float(
                np.linalg.norm(columns[:, i] / norms[i] - columns[:, j] / norms[j])
            )
            if distance <= response_tolerance:
                pairs.append(
                    {
                        "indices": [i, j],
                        "basis_type": family,
                        "response_distance": distance,
                    }
                )
    return pairs


def reduce_y_admittance_urn(
    model: YAdmittanceURN,
    freqs_hz,
    z_data,
    *,
    uncertainty_ohm,
    refit_epochs=6000,
    refit_restarts=3,
    response_tolerance=0.01,
):
    """Refit a nested path down to one basis; return smallest acceptable trial.

    ``uncertainty_ohm`` is a positive scalar or sample-aligned complex-magnitude
    error budget. Acceptance requires every |Zfit-Zdata|/budget <= 1. It is not
    an S-domain tolerance or a probabilistic confidence interval. The caller
    must supply its measurement meaning; no precision is inferred from data.
    Duplicate candidates are removed first, otherwise the least influential
    basis is removed. Rejected trials remain in the trace and are never labeled
    acceptable. This checks training samples, not independent generalization.
    The input model is not mutated; trace contains no raw measurement samples.
    """
    freqs = np.asarray(freqs_hz, dtype=float)
    target = np.asarray(z_data, dtype=complex)
    if (
        freqs.ndim != 1
        or len(freqs) < 2
        or target.shape != freqs.shape
        or not np.all(np.isfinite(freqs) & (freqs > 0))
        or len(np.unique(freqs)) != len(freqs)
        or not np.all(np.isfinite(target))
        or np.any(np.abs(target) == 0)
    ):
        raise ValueError(
            "need aligned finite nonzero impedance and unique positive frequencies"
        )
    budget = np.asarray(uncertainty_ohm, dtype=float)
    if budget.ndim != 0 and budget.shape != target.shape:
        raise ValueError("uncertainty_ohm must be scalar or sample-aligned")
    if not np.all(np.isfinite(budget) & (budget > 0)):
        raise ValueError("uncertainty_ohm must be finite and positive")
    for value in (refit_epochs, refit_restarts):
        if isinstance(value, bool) or not isinstance(value, int) or value < 1:
            raise ValueError("refit epochs/restarts must be positive integers")
    current, selected = model, None
    trials = []
    removed = None
    while True:
        prediction = current.predict(freqs)
        if not np.all(np.isfinite(prediction)):
            raise RuntimeError("nonfinite reduction prediction")
        ratio = float(np.max(np.abs(prediction - target) / budget))
        pairs = redundant_y_bases(current, freqs, response_tolerance=response_tolerance)
        accepted = ratio <= 1.0
        trials.append(
            {
                "basis_count": current.config.total_basis_functions,
                "max_uncertainty_ratio": ratio,
                "accepted": accepted,
                "s_domain_rmse": s_domain_rmse(prediction, target, z0=model.z0),
                "log_component_rmse": log_component_rmse(
                    prediction,
                    target,
                    floor=model.z0 * model.config.log_loss_floor_relative,
                ),
                "redundant_pairs": pairs,
                "removed_basis": removed,
            }
        )
        if accepted:
            selected = current
        if current.config.total_basis_functions == 1:
            break
        active = current.active_bases(freqs, threshold=0.0)
        duplicate_indices = {i for pair in pairs for i in pair["indices"]}
        pool = [b for b in active if b.basis_index in duplicate_indices] or active
        victim = min(pool, key=lambda b: (b.importance, b.basis_index))
        removed = {
            "index": victim.basis_index,
            "basis_type": victim.basis_type,
            "reason": "redundant_response" if duplicate_indices else "output_ablation",
        }
        retained = [b for b in active if b.basis_index != victim.basis_index]
        config = replace(
            current.config,
            n_epochs=refit_epochs,
            n_restarts=refit_restarts,
            sparsity_weight=0.0,
            z0=current.z0,
            y_ref=current.y_ref,
            omega_ref=current.omega_ref,
        )
        current = refit_y_admittance_active_bases(
            freqs, target, retained, config, verbose=False
        )
    return selected, {
        "schema": "radia.yurn-reduction.v1",
        "status": "accepted" if selected is not None else "no_acceptable_candidate",
        "source_basis_count": model.config.total_basis_functions,
        "selected_basis_count": None
        if selected is None
        else selected.config.total_basis_functions,
        "selection": "smallest tested candidate within caller-supplied per-sample ohmic budget",
        "uncertainty_ohm_range": [float(np.min(budget)), float(np.max(budget))],
        "response_tolerance": response_tolerance,
        "refit_epochs": refit_epochs,
        "refit_restarts": refit_restarts,
        "independent_validation": False,
        "trials": trials,
    }
