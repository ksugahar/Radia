"""Solver-neutral frequency-response gates for symmetric conductors."""

from __future__ import annotations

import math
from collections.abc import Iterable, Sequence
from typing import Any, Mapping


def twin_conductor_skin_effect_frequency_gate(
    frequencies_hz: Iterable[float],
    resistance_ohm: Iterable[Sequence[float]],
    inductance_h: Iterable[Sequence[float]],
    *,
    symmetry_rtol: float = 5.0e-4,
) -> dict[str, object]:
    """Gate passive, symmetric R/L frequency trends and derived impedance."""

    frequencies = [float(value) for value in frequencies_hz]
    resistance = [[float(value) for value in row] for row in resistance_ohm]
    inductance = [[float(value) for value in row] for row in inductance_h]
    if len(frequencies) < 5 or len(resistance) != len(frequencies) or len(inductance) != len(frequencies):
        raise ValueError("frequency, resistance, and inductance rows must have the same length >= 5")
    if any(len(row) != 2 for row in resistance + inductance):
        raise ValueError("each resistance and inductance row must contain two conductors")
    if not math.isfinite(float(symmetry_rtol)) or symmetry_rtol < 0.0:
        raise ValueError("symmetry_rtol must be finite and non-negative")
    scalars = frequencies + [value for row in resistance + inductance for value in row]
    if not all(math.isfinite(value) for value in scalars):
        raise ValueError("frequency-response values must be finite")

    symmetry_errors = [
        abs(a - b) / max(abs(a), abs(b), 1.0e-300)
        for row in resistance + inductance
        for a, b in [row]
    ]
    impedance = [
        [math.hypot(r, 2.0 * math.pi * frequency * l) for r, l in zip(r_row, l_row)]
        for frequency, r_row, l_row in zip(frequencies, resistance, inductance)
    ]
    checks = {
        "frequency_strictly_increasing_positive": all(value > 0.0 for value in frequencies)
        and all(a < b for a, b in zip(frequencies, frequencies[1:])),
        "resistance_positive": all(value > 0.0 for row in resistance for value in row),
        "inductance_positive": all(value > 0.0 for row in inductance for value in row),
        "resistance_non_decreasing": all(
            all(a <= b for a, b in zip(series, series[1:]))
            for series in zip(*resistance)
        ),
        "inductance_non_increasing": all(
            all(a >= b for a, b in zip(series, series[1:]))
            for series in zip(*inductance)
        ),
        "twin_conductor_symmetry": max(symmetry_errors) <= float(symmetry_rtol),
        "skin_effect_resistance_growth": all(resistance[-1][index] >= 2.0 * resistance[0][index] for index in range(2)),
        "impedance_magnitude_increasing": all(
            all(a < b for a, b in zip(series, series[1:]))
            for series in zip(*impedance)
        ),
    }
    return {
        "policy": "twin_conductor_skin_effect_frequency_gate_v1",
        "status": "ok" if all(checks.values()) else "needs_attention",
        "checks": checks,
        "issues": [name for name, ok in checks.items() if not ok],
        "metrics": {
            "frequency_count": len(frequencies),
            "maximum_twin_relative_mismatch": max(symmetry_errors),
            "resistance_growth_ratio": [resistance[-1][index] / resistance[0][index] for index in range(2)],
            "inductance_retention_ratio": [inductance[-1][index] / inductance[0][index] for index in range(2)],
            "impedance_magnitude_ohm": impedance,
        },
        "lesson": (
            "A passive conductor frequency sweep should keep R and L positive, show non-decreasing R "
            "and non-increasing L as skin/proximity effects develop, preserve geometric twin symmetry, "
            "and yield increasing |R+j omega L|."
        ),
    }


def homogenized_bundle_impedance_comparison_gate(
    rows: Iterable[Mapping[str, Any]],
    *,
    resistance_rtol: float = 0.03,
    inductance_rtol: float = 0.005,
    impedance_rtol: float = 0.01,
    observable_rtol: float = 1.0e-10,
    minimum_element_reduction: float = 5.0,
    minimum_speedup: float = 5.0,
) -> dict[str, object]:
    """Compare a homogenized stranded bundle with an explicit reference.

    The gate verifies passive complex impedance, recomputes ``Z=V/I`` and
    ``L=Im(Z)/omega``, then balances approximation error against mesh and solve
    cost. It is solver-neutral and suitable for round-wire homogenization,
    litz-wire surrogates, and explicit-strand reference models.
    """

    items = [dict(row) for row in rows]
    if len(items) != 2:
        raise ValueError("rows must contain one homogenized and one explicit_reference model")
    by_role = {str(row.get("model_role") or "").strip().lower(): row for row in items}
    if set(by_role) != {"homogenized", "explicit_reference"}:
        raise ValueError("model_role values must be homogenized and explicit_reference")

    normalized: dict[str, dict[str, Any]] = {}
    for role, row in by_role.items():
        try:
            frequency = float(row["frequency_hz"])
            current_pair = [float(value) for value in row["current_a_complex"]]
            voltage_pair = [float(value) for value in row["voltage_v_complex"]]
            resistance = float(row["resistance_ohm"])
            inductance = float(row["inductance_h"])
            elements = int(row["element_count"])
            solve_time = float(row["solve_time_s"])
        except (KeyError, TypeError, ValueError) as exc:
            raise ValueError(f"{role} row has missing or invalid observables") from exc
        if len(current_pair) != 2 or len(voltage_pair) != 2:
            raise ValueError("complex current and voltage must be [real, imag]")
        scalars = [frequency, *current_pair, *voltage_pair, resistance, inductance, solve_time]
        if not all(math.isfinite(value) for value in scalars):
            raise ValueError("all observables must be finite")
        current = complex(*current_pair)
        voltage = complex(*voltage_pair)
        if frequency <= 0.0 or abs(current) == 0.0 or elements <= 0 or solve_time <= 0.0:
            raise ValueError("frequency, current magnitude, element count, and solve time must be positive")
        impedance = voltage / current
        derived_l = impedance.imag / (2.0 * math.pi * frequency)
        normalized[role] = {
            "frequency_hz": frequency,
            "current": current,
            "impedance": impedance,
            "resistance_ohm": resistance,
            "inductance_h": inductance,
            "derived_resistance_ohm": impedance.real,
            "derived_inductance_h": derived_l,
            "element_count": elements,
            "solve_time_s": solve_time,
        }

    approximate = normalized["homogenized"]
    reference = normalized["explicit_reference"]
    resistance_error = abs(approximate["resistance_ohm"] - reference["resistance_ohm"]) / abs(reference["resistance_ohm"])
    inductance_error = abs(approximate["inductance_h"] - reference["inductance_h"]) / abs(reference["inductance_h"])
    impedance_error = abs(approximate["impedance"] - reference["impedance"]) / abs(reference["impedance"])
    element_reduction = reference["element_count"] / approximate["element_count"]
    speedup = reference["solve_time_s"] / approximate["solve_time_s"]

    observable_errors = {}
    for role, row in normalized.items():
        observable_errors[role] = {
            "resistance_relative": abs(row["derived_resistance_ohm"] - row["resistance_ohm"]) / max(abs(row["resistance_ohm"]), 1.0e-300),
            "inductance_relative": abs(row["derived_inductance_h"] - row["inductance_h"]) / max(abs(row["inductance_h"]), 1.0e-300),
        }

    checks = {
        "frequency_matches": approximate["frequency_hz"] == reference["frequency_hz"],
        "current_phasor_matches": approximate["current"] == reference["current"],
        "passive_positive_resistance": all(row["resistance_ohm"] > 0.0 for row in normalized.values()),
        "positive_series_inductance": all(row["inductance_h"] > 0.0 for row in normalized.values()),
        "reported_observables_match_voltage_current": all(
            error <= float(observable_rtol)
            for errors in observable_errors.values()
            for error in errors.values()
        ),
        "homogenized_resistance_accurate": resistance_error <= float(resistance_rtol),
        "homogenized_inductance_accurate": inductance_error <= float(inductance_rtol),
        "homogenized_complex_impedance_accurate": impedance_error <= float(impedance_rtol),
        "explicit_reference_has_more_elements": element_reduction >= float(minimum_element_reduction),
        "homogenized_model_is_faster": speedup >= float(minimum_speedup),
    }
    return {
        "policy": "homogenized_bundle_impedance_comparison_gate_v1",
        "status": "ok" if all(checks.values()) else "needs_attention",
        "checks": checks,
        "issues": [name for name, ok in checks.items() if not ok],
        "metrics": {
            "resistance_relative_error": resistance_error,
            "inductance_relative_error": inductance_error,
            "complex_impedance_relative_error": impedance_error,
            "element_count_reduction": element_reduction,
            "solve_time_speedup": speedup,
            "observable_reconstruction_errors": observable_errors,
        },
        "lesson": (
            "A strand homogenization is useful only when passive complex impedance agrees with an explicit "
            "reference and the saved element/time reduction is demonstrated. Compare R, L, and complex Z; "
            "a small |Z| error can otherwise hide a material resistance error."
        ),
    }


def opposed_busbar_skin_force_gate(
    rows: Iterable[Mapping[str, Any]],
    *,
    conductor_thickness_mm: float,
    conductivity_s_per_m: float,
    commanded_current_a: float,
    replay_rtol: float = 1.0e-12,
    identity_rtol: float = 5.0e-8,
    force_balance_rtol: float = 5.0e-5,
) -> dict[str, object]:
    """Gate skin/proximity, phasor power, and action-reaction force together."""

    items = [dict(row) for row in rows]
    thickness = float(conductor_thickness_mm)
    conductivity = float(conductivity_s_per_m)
    current = float(commanded_current_a)
    limits = [float(replay_rtol), float(identity_rtol), float(force_balance_rtol)]
    if len(items) != 8:
        raise ValueError("rows must contain two four-frequency replays")
    if thickness <= 0.0 or conductivity <= 0.0 or current <= 0.0:
        raise ValueError("thickness, conductivity, and commanded current must be positive")
    if not all(math.isfinite(value) and value >= 0.0 for value in limits):
        raise ValueError("relative tolerances must be finite and nonnegative")

    numeric_keys = (
        "frequency_hz",
        "skin_depth_mm",
        "ac_resistance_ohm",
        "circuit_current_re_a",
        "circuit_current_im_a",
        "voltage_re_v",
        "voltage_im_v",
        "flux_re_wb_turn",
        "total_loss_w",
        "top_loss_w",
        "bottom_loss_w",
        "top_current_re_a",
        "top_current_im_a",
        "bottom_current_re_a",
        "bottom_current_im_a",
        "top_lorentz_y_re_n",
        "bottom_lorentz_y_re_n",
        "energy_j",
        "coenergy_j",
        "inner_j_magnitude_ma_m2",
        "center_j_magnitude_ma_m2",
        "outer_j_magnitude_ma_m2",
    )
    normalized = []
    for index, row in enumerate(items):
        try:
            parsed = {key: float(row[key]) for key in numeric_keys}
            parsed["replay"] = int(row["replay"])
            parsed["node_count"] = int(row["node_count"])
            parsed["element_count"] = int(row["element_count"])
        except (KeyError, TypeError, ValueError) as exc:
            raise ValueError(f"row {index} has missing or invalid observables") from exc
        if not all(math.isfinite(value) for key, value in parsed.items() if key not in {"replay", "node_count", "element_count"}):
            raise ValueError(f"row {index} contains non-finite observables")
        normalized.append(parsed)

    by_replay = {
        replay: sorted(
            [row for row in normalized if row["replay"] == replay],
            key=lambda row: row["frequency_hz"],
        )
        for replay in (1, 2)
    }
    first, second = by_replay[1], by_replay[2]
    observable_keys = [*numeric_keys, "node_count", "element_count"]

    def relative(left: float, right: float) -> float:
        return abs(left - right) / max(abs(left), abs(right), 1.0e-300)

    replay_errors = [
        relative(left[key], right[key])
        for left, right in zip(first, second)
        for key in observable_keys
    ]
    skin_errors = []
    power_errors = []
    faraday_errors = []
    current_errors = []
    force_errors = []
    loss_split_errors = []
    loss_symmetry_errors = []
    energy_errors = []
    for row in normalized:
        omega = 2.0 * math.pi * row["frequency_hz"]
        analytic_skin_mm = math.sqrt(2.0 / (omega * 4.0e-7 * math.pi * conductivity)) * 1000.0
        skin_errors.append(relative(row["skin_depth_mm"], analytic_skin_mm))
        power_errors.append(relative(row["total_loss_w"], 0.5 * row["voltage_re_v"] * row["circuit_current_re_a"]))
        faraday_errors.append(relative(row["voltage_im_v"], omega * row["flux_re_wb_turn"]))
        top_current = complex(row["top_current_re_a"], row["top_current_im_a"])
        bottom_current = complex(row["bottom_current_re_a"], row["bottom_current_im_a"])
        current_errors.append(abs(top_current + bottom_current) / max(abs(top_current), abs(bottom_current), 1.0e-300))
        force_errors.append(abs(row["top_lorentz_y_re_n"] + row["bottom_lorentz_y_re_n"]) / max(abs(row["top_lorentz_y_re_n"]), abs(row["bottom_lorentz_y_re_n"]), 1.0e-300))
        loss_split_errors.append(relative(row["total_loss_w"], row["top_loss_w"] + row["bottom_loss_w"]))
        loss_symmetry_errors.append(relative(row["top_loss_w"], row["bottom_loss_w"]))
        energy_errors.append(relative(row["energy_j"], row["coenergy_j"]))

    frequencies = [row["frequency_hz"] for row in first]
    resistance = [row["ac_resistance_ohm"] for row in first]
    low, high = first[0], first[-1]
    checks = {
        "two_complete_replays": len(first) == len(second) == 4,
        "shared_strictly_increasing_frequency_axis": [row["frequency_hz"] for row in second] == frequencies and all(left < right for left, right in zip(frequencies, frequencies[1:])),
        "replay_observables_match": max(replay_errors) <= float(replay_rtol),
        "positive_stable_mesh_inventory": len({(row["node_count"], row["element_count"]) for row in normalized}) == 1 and min(row["node_count"] for row in normalized) > 0 and min(row["element_count"] for row in normalized) > 0,
        "commanded_circuit_current_preserved": max(abs(row["circuit_current_re_a"] - current) for row in normalized) <= 1.0e-12 and max(abs(row["circuit_current_im_a"]) for row in normalized) <= 1.0e-12,
        "opposed_conductor_current_closure": max(current_errors) <= 5.0e-7,
        "analytic_skin_depth_reproduced": max(skin_errors) <= float(replay_rtol),
        "phasor_real_power_identity": max(power_errors) <= float(identity_rtol),
        "faraday_voltage_flux_identity": max(faraday_errors) <= float(identity_rtol),
        "positive_loss_and_exact_partition": min(row["total_loss_w"] for row in normalized) > 0.0 and max(loss_split_errors) <= 1.0e-12,
        "symmetric_conductor_loss": max(loss_symmetry_errors) <= 1.0e-5,
        "linear_energy_coenergy_identity": max(energy_errors) <= float(replay_rtol),
        "lorentz_action_reaction_closure": all(row["top_lorentz_y_re_n"] > 0.0 and row["bottom_lorentz_y_re_n"] < 0.0 for row in normalized) and max(force_errors) <= float(force_balance_rtol),
        "ac_resistance_monotonic_and_doubled": all(left < right for left, right in zip(resistance, resistance[1:])) and resistance[-1] / resistance[0] >= 2.0,
        "low_frequency_density_nearly_uniform": low["inner_j_magnitude_ma_m2"] / low["center_j_magnitude_ma_m2"] <= 1.01 and low["inner_j_magnitude_ma_m2"] / low["outer_j_magnitude_ma_m2"] <= 1.01,
        "high_frequency_inner_face_crowding": high["skin_depth_mm"] < thickness and high["inner_j_magnitude_ma_m2"] / high["center_j_magnitude_ma_m2"] >= 3.0 and high["inner_j_magnitude_ma_m2"] / high["outer_j_magnitude_ma_m2"] >= 5.0,
    }
    issues = [name for name, ok in checks.items() if not ok]
    return {
        "policy": "opposed_busbar_skin_force_gate_v1",
        "status": "ok" if not issues else "needs_attention",
        "checks": checks,
        "issues": issues,
        "metrics": {
            "maximum_replay_relative_error": max(replay_errors),
            "maximum_skin_depth_relative_error": max(skin_errors),
            "maximum_power_identity_relative_error": max(power_errors),
            "maximum_faraday_identity_relative_error": max(faraday_errors),
            "maximum_current_closure_relative_error": max(current_errors),
            "maximum_force_closure_relative_error": max(force_errors),
            "resistance_growth_ratio": resistance[-1] / resistance[0],
            "inner_to_center_density_high_frequency": high["inner_j_magnitude_ma_m2"] / high["center_j_magnitude_ma_m2"],
            "inner_to_outer_density_high_frequency": high["inner_j_magnitude_ma_m2"] / high["outer_j_magnitude_ma_m2"],
        },
        "lesson": (
            "An opposed-conductor AC sweep should close current and Lorentz action-reaction, satisfy phasor power and Faraday identities, reproduce analytic skin depth, and resolve inner-face proximity crowding while resistance rises."
        ),
    }
