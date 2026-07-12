"""Solver-neutral validation gate for reduced proximity-effect bundle models."""

from __future__ import annotations

import json
import math
from typing import Any


def _relative_error(value: float, reference: float) -> float:
    return abs(value - reference) / max(abs(reference), 1.0e-30)


def litz_proximity_approximation_pair_gate(summary_json: str) -> dict[str, Any]:
    """Check a reduced conductor-bundle model against an explicit reference.

    The contract is deliberately solver-neutral. Both models must report the
    same prescribed current and ampere-turns, complex impedance, total loss,
    complex-power closure, field probes, mesh size, and solve time.
    """

    policy = "litz_proximity_approximation_pair_gate_v1"
    try:
        payload = json.loads(summary_json)
    except (json.JSONDecodeError, TypeError) as exc:
        return {"policy": policy, "status": "invalid_input", "error": str(exc)}
    if not isinstance(payload, dict):
        return {"policy": policy, "status": "invalid_input", "error": "summary must be an object"}

    models = payload.get("models")
    if not isinstance(models, dict) or set(models) != {"approximate", "exact"}:
        return {
            "policy": policy,
            "status": "needs_attention",
            "errors": ["models must contain exactly approximate and exact"],
        }

    parsed: dict[str, dict[str, Any]] = {}
    errors: list[str] = []
    for name in ("approximate", "exact"):
        row = models[name]
        try:
            impedance = complex(float(row["impedance_ohm"]["real"]), float(row["impedance_ohm"]["imag"]))
            probes = {
                str(item["id"]): float(item["b_abs_t"])
                for item in row["field_samples"]
            }
            values = {
                "current_a": float(row["current_a"]),
                "total_current_a_turn": float(row["total_current_a_turn"]),
                "loss_w": float(row["total_loss_w"]),
                "power_w": float(row["half_real_vi_star_w"]),
                "elements": int(row["element_count"]),
                "solve_s": float(row["solve_time_s"]),
            }
        except (KeyError, TypeError, ValueError) as exc:
            errors.append(f"{name} observables invalid: {exc}")
            continue
        finite = [impedance.real, impedance.imag, *values.values(), *probes.values()]
        if not all(math.isfinite(value) for value in finite):
            errors.append(f"{name} observables must be finite")
        parsed[name] = {"impedance": impedance, "probes": probes, **values}

    if len(parsed) != 2:
        return {"policy": policy, "status": "needs_attention", "errors": errors}

    approx = parsed["approximate"]
    exact = parsed["exact"]
    common_probes = set(approx["probes"]) & set(exact["probes"])
    field_errors = [
        _relative_error(approx["probes"][probe], exact["probes"][probe])
        for probe in sorted(common_probes)
    ]
    metrics = {
        "total_loss_relative_error": _relative_error(approx["loss_w"], exact["loss_w"]),
        "resistance_relative_error": _relative_error(approx["impedance"].real, exact["impedance"].real),
        "reactance_relative_error": _relative_error(approx["impedance"].imag, exact["impedance"].imag),
        "complex_impedance_relative_error": abs(approx["impedance"] - exact["impedance"])
        / max(abs(exact["impedance"]), 1.0e-30),
        "max_field_sample_relative_error": max(field_errors, default=math.inf),
        "element_count_ratio_exact_over_approx": exact["elements"] / max(approx["elements"], 1),
        "solve_time_ratio_exact_over_approx": exact["solve_s"] / max(approx["solve_s"], 1.0e-30),
        "approximate_power_identity_relative_error": _relative_error(approx["loss_w"], approx["power_w"]),
        "exact_power_identity_relative_error": _relative_error(exact["loss_w"], exact["power_w"]),
    }
    checks = {
        "schema_recorded": payload.get("schema") == "litz.proximity-approximation-pair.v1",
        "frequency_recorded": float(payload.get("frequency_hz", 0.0)) > 0.0,
        "same_one_amp_excitation": all(abs(row["current_a"] - 1.0) <= 1.0e-9 for row in parsed.values()),
        "same_114_ampere_turns": all(abs(row["total_current_a_turn"] - 114.0) <= 1.0e-6 for row in parsed.values()),
        "passive_positive_impedance_and_loss": all(
            row["impedance"].real > 0.0 and row["impedance"].imag > 0.0 and row["loss_w"] > 0.0
            for row in parsed.values()
        ),
        "three_matching_field_probes": len(common_probes) == 3
        and len(approx["probes"]) == 3
        and len(exact["probes"]) == 3,
        "complex_power_closes_total_loss": metrics["approximate_power_identity_relative_error"] <= 2.0e-4
        and metrics["exact_power_identity_relative_error"] <= 2.0e-4,
        "loss_approximation_within_3pct": metrics["total_loss_relative_error"] <= 0.03,
        "resistance_approximation_within_3pct": metrics["resistance_relative_error"] <= 0.03,
        "reactance_approximation_within_half_percent": metrics["reactance_relative_error"] <= 0.005,
        "complex_impedance_within_1pct": metrics["complex_impedance_relative_error"] <= 0.01,
        "sampled_field_within_4pct": metrics["max_field_sample_relative_error"] <= 0.04,
        "explicit_model_is_at_least_5x_larger": metrics["element_count_ratio_exact_over_approx"] >= 5.0,
        "reduced_model_is_at_least_5x_faster": metrics["solve_time_ratio_exact_over_approx"] >= 5.0,
    }
    failed = [name for name, ok in checks.items() if not ok]
    return {
        "policy": policy,
        "status": "ok" if not errors and not failed else "needs_attention",
        "checks": checks,
        "metrics": metrics,
        "errors": errors + failed,
        "lesson": (
            "Accept a reduced high-frequency conductor bundle only after total loss, complex power, "
            "complex impedance, and spatial field probes agree with an explicit conductor reference; "
            "a speedup alone is not validation."
        ),
    }
