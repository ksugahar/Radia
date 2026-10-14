"""Solver-neutral gate for two-terminal stationary-current result packages."""

from __future__ import annotations

import math
from collections.abc import Mapping, Sequence
from typing import Any


def _finite(value: object, name: str) -> float:
    try:
        parsed = float(value)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"{name} must be numeric") from exc
    if not math.isfinite(parsed):
        raise ValueError(f"{name} must be finite")
    return parsed


def _rows(value: object, name: str, minimum: int) -> list[Mapping[str, object]]:
    if not isinstance(value, Sequence) or isinstance(value, (str, bytes)):
        raise ValueError(f"{name} must be an array")
    rows = list(value)
    if len(rows) < minimum or any(not isinstance(row, Mapping) for row in rows):
        raise ValueError(f"{name} must contain at least {minimum} objects")
    return rows


def _pair(value: object, name: str) -> tuple[float, float]:
    if not isinstance(value, Sequence) or isinstance(value, (str, bytes)) or len(value) != 2:
        raise ValueError(f"{name} must contain exactly two values")
    return _finite(value[0], f"{name}[0]"), _finite(value[1], f"{name}[1]")


def _relative_error(left: float, right: float) -> float:
    return abs(left - right) / max(abs(left), abs(right), 1.0e-300)


def two_terminal_dc_conduction_power_gate(
    summary: Mapping[str, object],
) -> dict[str, Any]:
    """Recompute current, power, loss, adaptive, and fresh-run closure."""

    if not isinstance(summary, Mapping):
        raise ValueError("summary must be an object")
    runs = _rows(summary.get("runs"), "runs", 2)
    tolerances = summary.get("tolerances", {})
    if not isinstance(tolerances, Mapping):
        raise ValueError("tolerances must be an object")
    limits = {
        "current_balance": _finite(
            tolerances.get("maximum_relative_current_balance_error", 1.0e-5),
            "maximum_relative_current_balance_error",
        ),
        "power_loss": _finite(
            tolerances.get("maximum_relative_terminal_power_loss_error", 1.0e-5),
            "maximum_relative_terminal_power_loss_error",
        ),
        "loss_parts": _finite(
            tolerances.get("maximum_relative_loss_decomposition_error", 1.0e-8),
            "maximum_relative_loss_decomposition_error",
        ),
        "adaptive_error": _finite(
            tolerances.get("maximum_final_adaptive_relative_error", 1.0e-3),
            "maximum_final_adaptive_relative_error",
        ),
        "adaptive_change": _finite(
            tolerances.get("maximum_final_adaptive_power_change", 1.0e-3),
            "maximum_final_adaptive_power_change",
        ),
        "repeat": _finite(
            tolerances.get("maximum_fresh_run_relative_difference", 1.0e-9),
            "maximum_fresh_run_relative_difference",
        ),
    }
    if any(value < 0.0 for value in limits.values()):
        raise ValueError("tolerances must be nonnegative")

    metrics = []
    repeat_vectors: list[list[float]] = []
    run_ids: list[str] = []
    for index, run in enumerate(runs):
        run_id = str(run.get("run_id") or "").strip()
        currents = _pair(run.get("port_currents_a"), f"runs[{index}].port_currents_a")
        potentials = _pair(run.get("port_potentials_v"), f"runs[{index}].port_potentials_v")
        integrated_loss = _finite(
            run.get("integrated_loss_w"), f"runs[{index}].integrated_loss_w"
        )
        loss_parts_value = run.get("loss_parts_w")
        if not isinstance(loss_parts_value, Sequence) or isinstance(
            loss_parts_value, (str, bytes)
        ) or not loss_parts_value:
            raise ValueError(f"runs[{index}].loss_parts_w must be a nonempty array")
        loss_parts = [
            _finite(value, f"runs[{index}].loss_parts_w[{part}]")
            for part, value in enumerate(loss_parts_value)
        ]
        adaptive = _rows(run.get("adaptive_rows"), f"runs[{index}].adaptive_rows", 2)
        parsed_adaptive = []
        for row_index, row in enumerate(adaptive):
            parsed_adaptive.append(
                {
                    "mesh_cells": int(
                        _finite(row.get("mesh_cells"), f"adaptive[{row_index}].mesh_cells")
                    ),
                    "degrees_of_freedom": int(
                        _finite(
                            row.get("degrees_of_freedom"),
                            f"adaptive[{row_index}].degrees_of_freedom",
                        )
                    ),
                    "power_w": _finite(row.get("power_w"), f"adaptive[{row_index}].power_w"),
                }
            )
        final_error = _finite(
            run.get("final_adaptive_relative_error"),
            f"runs[{index}].final_adaptive_relative_error",
        )
        through_current = 0.5 * (abs(currents[0]) + abs(currents[1]))
        terminal_voltage = abs(potentials[1] - potentials[0])
        terminal_power = sum(v * i for v, i in zip(potentials, currents))
        resistance = terminal_voltage / max(through_current, 1.0e-300)
        current_balance = abs(sum(currents)) / max(through_current, 1.0e-300)
        power_loss_error = _relative_error(terminal_power, integrated_loss)
        loss_parts_error = _relative_error(sum(loss_parts), integrated_loss)
        final_power_error = _relative_error(parsed_adaptive[-1]["power_w"], integrated_loss)
        final_power_change = _relative_error(
            parsed_adaptive[-1]["power_w"], parsed_adaptive[-2]["power_w"]
        )
        metric = {
            "run_id": run_id,
            "solver_complete": run.get("solver_complete") is True,
            "through_current_a": through_current,
            "terminal_voltage_v": terminal_voltage,
            "terminal_power_w": terminal_power,
            "integrated_loss_w": integrated_loss,
            "effective_resistance_ohm": resistance,
            "relative_current_balance_error": current_balance,
            "relative_terminal_power_loss_error": power_loss_error,
            "relative_loss_decomposition_error": loss_parts_error,
            "adaptive_pass_count": len(parsed_adaptive),
            "final_mesh_cells": parsed_adaptive[-1]["mesh_cells"],
            "final_degrees_of_freedom": parsed_adaptive[-1]["degrees_of_freedom"],
            "final_adaptive_relative_error": final_error,
            "final_adaptive_power_relative_change": final_power_change,
            "final_adaptive_power_loss_error": final_power_error,
            "mesh_cells_strictly_increase": all(
                left["mesh_cells"] < right["mesh_cells"]
                for left, right in zip(parsed_adaptive, parsed_adaptive[1:])
            ),
            "degrees_of_freedom_strictly_increase": all(
                left["degrees_of_freedom"] < right["degrees_of_freedom"]
                for left, right in zip(parsed_adaptive, parsed_adaptive[1:])
            ),
        }
        metrics.append(metric)
        run_ids.append(run_id)
        repeat_vectors.append(
            [
                through_current,
                terminal_voltage,
                terminal_power,
                integrated_loss,
                resistance,
                float(metric["final_mesh_cells"]),
                float(metric["final_degrees_of_freedom"]),
                final_error,
            ]
        )

    repeat_errors = [
        max(_relative_error(left, right) for left, right in zip(repeat_vectors[0], vector))
        for vector in repeat_vectors[1:]
    ]
    maximum_repeat_error = max(repeat_errors, default=math.inf)
    checks = {
        "two_or_more_distinct_complete_runs": len(set(run_ids)) == len(runs)
        and all(run_ids)
        and all(row["solver_complete"] for row in metrics),
        "opposed_terminal_currents_close": max(
            row["relative_current_balance_error"] for row in metrics
        )
        <= limits["current_balance"],
        "terminal_power_matches_integrated_joule_loss": max(
            row["relative_terminal_power_loss_error"] for row in metrics
        )
        <= limits["power_loss"],
        "loss_decomposition_closes": max(
            row["relative_loss_decomposition_error"] for row in metrics
        )
        <= limits["loss_parts"],
        "effective_resistance_is_finite_positive": all(
            math.isfinite(row["effective_resistance_ohm"])
            and row["effective_resistance_ohm"] > 0.0
            for row in metrics
        ),
        "adaptive_mesh_and_dof_strictly_increase": all(
            row["mesh_cells_strictly_increase"]
            and row["degrees_of_freedom_strictly_increase"]
            for row in metrics
        ),
        "final_adaptive_error_is_bounded": max(
            row["final_adaptive_relative_error"] for row in metrics
        )
        <= limits["adaptive_error"],
        "adaptive_power_converges_to_integrated_loss": max(
            max(
                row["final_adaptive_power_relative_change"],
                row["final_adaptive_power_loss_error"],
            )
            for row in metrics
        )
        <= limits["adaptive_change"],
        "fresh_runs_repeat_scalar_observables": maximum_repeat_error <= limits["repeat"],
    }
    accepted = all(checks.values())
    return {
        "policy": "two_terminal_dc_conduction_power_gate_v1",
        "status": "ok" if accepted else "needs_attention",
        "solver_ready": accepted,
        "checks": checks,
        "issues": [name for name, ok in checks.items() if not ok],
        "limits": limits,
        "metrics": {
            "run_count": len(runs),
            "per_run": metrics,
            "maximum_fresh_run_relative_difference": maximum_repeat_error,
        },
        "lesson": (
            "A stationary-current result is credible only when opposed terminal currents close, "
            "terminal VI power equals integrated Joule loss, the loss partition closes, adaptive "
            "power converges, and a fresh solve reproduces the scalar observables."
        ),
    }
