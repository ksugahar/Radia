"""Validation gate for alternate volume and surface eddy-loss formulations."""
from __future__ import annotations

import math
from collections.abc import Mapping
from typing import Any


def _finite_positive(row: Mapping[str, Any], name: str) -> float:
    try:
        value = float(row[name])
    except (KeyError, TypeError, ValueError) as exc:
        raise ValueError(f"{name} must be a finite positive number") from exc
    if not math.isfinite(value) or value <= 0.0:
        raise ValueError(f"{name} must be a finite positive number")
    return value


def alternate_eddy_loss_formulation_gate(
    summary: Mapping[str, Any],
    *,
    identity_rtol: float = 1.0e-10,
    rerun_rtol: float = 1.0e-6,
) -> dict[str, object]:
    """Gate resolved-volume and surface-impedance results as alternatives.

    The two formulations may use different meshes, selections, and solution
    datasets. Their losses are therefore comparison observables, not additive
    contributions to one total.
    """
    if not isinstance(summary, Mapping):
        raise ValueError("summary must be a mapping")
    volume = summary.get("volume")
    surface = summary.get("surface")
    if not isinstance(volume, Mapping) or not isinstance(surface, Mapping):
        raise ValueError("volume and surface must be mappings")
    frequency = _finite_positive(summary, "frequency_hz")
    volume_native = _finite_positive(volume, "native_loss_w")
    volume_builtin = _finite_positive(volume, "builtin_integral_w")
    volume_jdot_e = _finite_positive(volume, "jdot_e_integral_w")
    surface_native = _finite_positive(surface, "native_loss_w")
    surface_builtin = _finite_positive(surface, "builtin_integral_w")
    try:
        rerun_change = float(surface["rerun_relative_change"])
        volume_dim = int(volume["selection_dimension"])
        surface_dim = int(surface["selection_dimension"])
    except (KeyError, TypeError, ValueError) as exc:
        raise ValueError("selection dimensions and rerun_relative_change are required") from exc
    if not math.isfinite(rerun_change) or rerun_change < 0.0:
        raise ValueError("rerun_relative_change must be finite and nonnegative")

    rel = lambda left, right: abs(left - right) / max(abs(right), 1.0e-300)
    volume_dataset = str(volume.get("dataset_id") or "")
    surface_dataset = str(surface.get("dataset_id") or "")
    volume_solution = str(volume.get("solution_id") or "")
    surface_solution = str(surface.get("solution_id") or "")
    checks = {
        "frequency_positive": frequency > 0.0,
        "dataset_ids_recorded_and_distinct": bool(volume_dataset and surface_dataset)
        and volume_dataset != surface_dataset,
        "solution_ids_recorded_and_distinct": bool(volume_solution and surface_solution)
        and volume_solution != surface_solution,
        "volume_selection_is_domain": volume_dim == 3,
        "surface_selection_is_boundary": surface_dim == 2,
        "volume_builtin_matches_native": rel(volume_builtin, volume_native) <= identity_rtol,
        "volume_jdot_e_identity": rel(volume_jdot_e, volume_builtin) <= identity_rtol,
        "surface_builtin_matches_native": rel(surface_builtin, surface_native) <= identity_rtol,
        "surface_rerun_reproducible": rerun_change <= rerun_rtol,
        "alternate_formulations_not_added": summary.get("combine_requested") is False,
    }
    return {
        "policy": "alternate_eddy_loss_formulation_gate_v1",
        "status": "ok" if all(checks.values()) else "needs_attention",
        "checks": checks,
        "issues": [name for name, ok in checks.items() if not ok],
        "metrics": {
            "frequency_hz": frequency,
            "surface_to_volume_loss_ratio": surface_native / volume_native,
            "volume_jdot_e_relative_error": rel(volume_jdot_e, volume_builtin),
            "surface_rerun_relative_change": rerun_change,
        },
        "lesson": (
            "Resolved-volume Joule loss and a surface-impedance loss on distinct "
            "solution datasets are alternate formulations. Validate each on its "
            "own geometric selection; compare them diagnostically and never add them."
        ),
    }
