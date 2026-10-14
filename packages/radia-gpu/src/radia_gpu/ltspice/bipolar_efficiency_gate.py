"""Solver-neutral dual-output converter efficiency gate."""

from __future__ import annotations

import math
from collections.abc import Mapping, Sequence


def bipolar_converter_efficiency_gate(summary: Mapping[str, object]) -> dict[str, object]:
    """Gate signed input power, dual-output balance, and efficiency closure."""

    units = summary.get("units")
    window = summary.get("measure_window_s")
    if not isinstance(units, Mapping):
        raise ValueError("units must be a mapping")
    if not isinstance(window, Sequence) or isinstance(window, (str, bytes)) or len(window) != 2:
        raise ValueError("measure_window_s must contain start and stop")

    def finite(name: str) -> float:
        try:
            value = float(summary[name])
        except (KeyError, TypeError, ValueError) as exc:
            raise ValueError(f"{name} must be finite") from exc
        if not math.isfinite(value):
            raise ValueError(f"{name} must be finite")
        return value

    start_s, stop_s = (float(value) for value in window)
    input_source_power_w = finite("input_source_power_w")
    output_1_power_w = finite("output_1_power_w")
    output_2_power_w = finite("output_2_power_w")
    reported_efficiency_percent = finite("reported_efficiency_percent")
    closure_limit = float(summary.get("max_efficiency_closure_percent", 1.0e-6))
    imbalance_limit = float(summary.get("max_output_imbalance_relative", 0.02))

    delivered_input_w = -input_source_power_w
    output_sum_w = output_1_power_w + output_2_power_w
    loss_w = delivered_input_w - output_sum_w
    recomputed_efficiency_percent = (
        100.0 * output_sum_w / delivered_input_w if delivered_input_w > 0.0 else math.nan
    )
    output_mean_w = output_sum_w / 2.0
    output_imbalance_relative = (
        abs(output_1_power_w - output_2_power_w) / output_mean_w
        if output_mean_w > 0.0
        else math.inf
    )
    efficiency_closure_percent = abs(
        reported_efficiency_percent - recomputed_efficiency_percent
    )

    checks = {
        "units_are_explicit": units.get("power") == "W"
        and units.get("efficiency") == "percent",
        "measure_window_is_ordered": math.isfinite(start_s)
        and math.isfinite(stop_s)
        and 0.0 <= start_s < stop_s,
        "source_branch_sign_means_delivered_power": input_source_power_w < 0.0,
        "both_outputs_deliver_positive_power": output_1_power_w > 0.0
        and output_2_power_w > 0.0,
        "power_balance_is_passive": delivered_input_w > 0.0
        and output_sum_w > 0.0
        and loss_w >= 0.0,
        "reported_efficiency_is_physical": 0.0 < reported_efficiency_percent <= 100.0,
        "reported_efficiency_recomputes": efficiency_closure_percent
        <= closure_limit,
        "dual_outputs_are_balanced": output_imbalance_relative <= imbalance_limit,
    }
    return {
        "policy": "bipolar_converter_efficiency_gate_v1",
        "status": "ok" if all(checks.values()) else "needs_attention",
        "checks": checks,
        "issues": [name for name, ok in checks.items() if not ok],
        "metrics": {
            "delivered_input_power_w": delivered_input_w,
            "output_power_sum_w": output_sum_w,
            "power_loss_w": loss_w,
            "recomputed_efficiency_percent": recomputed_efficiency_percent,
            "efficiency_closure_absolute_percent": efficiency_closure_percent,
            "output_power_imbalance_relative": output_imbalance_relative,
        },
        "lesson": (
            "Converter efficiency is credible only when the source-current sign is interpreted as delivered power, "
            "all output powers are positive, output power does not exceed input power, and the reported efficiency "
            "recomputes from the same late measurement window."
        ),
    }
