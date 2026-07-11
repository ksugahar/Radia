"""Solver-neutral gates for cyclic multi-terminal source sweeps."""

from __future__ import annotations

import math
from typing import Any


def cyclic_terminal_source_sweep_gate(
    summary: dict[str, Any],
    max_cyclic_relative_spread: float = 1.0e-8,
) -> dict[str, Any]:
    """Check cyclic symmetry without assuming two formulations are identical."""
    if not isinstance(summary, dict):
        raise ValueError("summary must be a mapping")
    formulations = summary.get("formulations") or []
    if not isinstance(formulations, list) or len(formulations) < 2:
        raise ValueError("at least two formulations are required")

    tolerance = float(max_cyclic_relative_spread)
    if not math.isfinite(tolerance) or tolerance < 0.0:
        raise ValueError("max_cyclic_relative_spread must be finite and nonnegative")

    checks: dict[str, bool] = {}
    metrics: list[dict[str, Any]] = []
    labels: list[str] = []
    for index, row in enumerate(formulations, start=1):
        if not isinstance(row, dict):
            raise ValueError(f"formulation {index} must be a mapping")
        label = str(row.get("label") or "").strip()
        charges = [float(value) for value in row.get("active_terminal_charge_c") or []]
        order = [int(value) for value in row.get("source_sweep_order") or []]
        terminal_count = len(charges)
        finite_positive = bool(charges) and all(
            math.isfinite(value) and value > 0.0 for value in charges
        )
        scale = max((abs(value) for value in charges), default=0.0)
        spread = (
            (max(charges) - min(charges)) / scale
            if finite_positive and scale > 0.0
            else math.inf
        )

        labels.append(label)
        prefix = f"formulation_{index}"
        checks[f"{prefix}_label_recorded"] = bool(label)
        checks[f"{prefix}_terminal_count_sufficient"] = terminal_count >= 3
        checks[f"{prefix}_order_is_permutation"] = (
            len(order) == terminal_count
            and sorted(order) == list(range(1, terminal_count + 1))
        )
        checks[f"{prefix}_charges_positive_finite"] = finite_positive
        checks[f"{prefix}_cyclic_symmetry"] = spread <= tolerance
        metrics.append(
            {
                "label": label,
                "terminal_count": terminal_count,
                "mean_active_charge_c": (
                    sum(charges) / terminal_count if finite_positive else None
                ),
                "cyclic_relative_spread": spread,
            }
        )

    checks["formulation_labels_distinct"] = (
        len(set(labels)) == len(labels) and all(labels)
    )
    means = [row["mean_active_charge_c"] for row in metrics]
    ratio = (
        max(means) / min(means)
        if means and all(value is not None and value > 0.0 for value in means)
        else math.inf
    )
    checks["cross_formulation_ratio_recorded_not_forced_to_unity"] = (
        math.isfinite(ratio) and ratio > 0.0
    )

    return {
        "policy": "cyclic_terminal_source_sweep_gate_v1",
        "status": "ok" if all(checks.values()) else "needs_attention",
        "checks": checks,
        "issues": [name for name, ok in checks.items() if not ok],
        "metrics": {
            "formulations": metrics,
            "maximum_to_minimum_formulation_ratio": ratio,
        },
        "tolerances": {"max_cyclic_relative_spread": tolerance},
        "notes": [
            "Record outer-solution order before extracting active-terminal charge.",
            "Cyclic symmetry is an internal invariant, not proof that two formulations share the same exterior model.",
            "Do not force absolute agreement when boundary and floating-potential contracts differ.",
        ],
    }
