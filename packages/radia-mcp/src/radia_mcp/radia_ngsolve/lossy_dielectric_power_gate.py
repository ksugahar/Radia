"""Solver-neutral complex-power gate for a lossy dielectric refinement study."""

from __future__ import annotations

import math
from collections.abc import Mapping, Sequence
from typing import Any


_DEFAULT_LIMITS = {
    "loss_ratio_relative": 1.0e-4,
    "reactive_energy_relative": 1.0e-4,
    "complex_real_relative": 5.0e-4,
    "complex_imag_relative": 1.0e-4,
    "complex_magnitude_relative": 2.0e-5,
    "observable_imag_relative": 2.0e-4,
    "last_pair_convergence_relative": 2.0e-3,
    "voltage_drop_span_relative": 1.0e-8,
}


def _finite(value: object, name: str, *, positive: bool = False) -> float:
    try:
        parsed = float(value)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"{name} must be numeric") from exc
    if not math.isfinite(parsed):
        raise ValueError(f"{name} must be finite")
    if positive and parsed <= 0.0:
        raise ValueError(f"{name} must be positive")
    return parsed


def _complex(value: object, name: str) -> tuple[complex, float]:
    if not isinstance(value, Mapping):
        raise ValueError(f"{name} must be an object with real, imag, and abs")
    real = _finite(value.get("real"), f"{name}.real")
    imag = _finite(value.get("imag"), f"{name}.imag")
    magnitude = _finite(value.get("abs"), f"{name}.abs")
    if magnitude < 0.0:
        raise ValueError(f"{name}.abs must be nonnegative")
    encoded = complex(real, imag)
    scale = max(abs(encoded), magnitude, 1.0e-300)
    if abs(abs(encoded) - magnitude) / scale > 1.0e-12:
        raise ValueError(f"{name}.abs is inconsistent with real and imag")
    return encoded, magnitude


def _relative(actual: float, reference: float) -> float:
    return abs(actual - reference) / max(abs(reference), 1.0e-300)


def lossy_dielectric_complex_power_refinement_gate(
    summary: Mapping[str, object],
) -> dict[str, Any]:
    """Recompute constitutive, energy, complex-power, and refinement checks."""

    if not isinstance(summary, Mapping):
        raise ValueError("summary must be an object")
    frequency = _finite(summary.get("frequency_Hz"), "frequency_Hz", positive=True)
    sigma = _finite(summary.get("sigma_S_per_m"), "sigma_S_per_m", positive=True)
    epsilon_r = _finite(summary.get("epsilon_r"), "epsilon_r", positive=True)
    epsilon_0 = _finite(
        summary.get("epsilon_0_F_per_m"), "epsilon_0_F_per_m", positive=True
    )
    rows_value = summary.get("rows")
    if not isinstance(rows_value, Sequence) or isinstance(rows_value, (str, bytes)):
        raise ValueError("rows must be an array")
    rows = list(rows_value)
    if len(rows) < 3 or any(not isinstance(row, Mapping) for row in rows):
        raise ValueError("rows must contain at least three result objects")

    overrides = summary.get("gate_tolerances", {})
    if not isinstance(overrides, Mapping):
        raise ValueError("gate_tolerances must be an object")
    limits = {
        name: _finite(overrides.get(name, default), name, positive=True)
        for name, default in _DEFAULT_LIMITS.items()
    }

    omega = 2.0 * math.pi * frequency
    expected_loss_ratio = sigma / (omega * epsilon_0 * epsilon_r)
    metrics: list[dict[str, float | int]] = []
    element_counts: list[int] = []
    mesh_sizes: list[float] = []
    powers: list[float] = []
    reactive_powers: list[float] = []
    voltage_drops: list[float] = []

    for index, row_value in enumerate(rows):
        row = row_value
        mesh_size = _finite(row.get("mesh_size_in"), f"rows[{index}].mesh_size_in", positive=True)
        element_count = int(_finite(row.get("element_count"), f"rows[{index}].element_count", positive=True))
        if element_count != _finite(row.get("element_count"), f"rows[{index}].element_count"):
            raise ValueError(f"rows[{index}].element_count must be an integer")
        p_complex, _ = _complex(row.get("real_power_W"), f"rows[{index}].real_power_W")
        q_complex, _ = _complex(row.get("reactive_power_var"), f"rows[{index}].reactive_power_var")
        s_complex, s_abs = _complex(row.get("apparent_power_VA"), f"rows[{index}].apparent_power_VA")
        w_complex, _ = _complex(
            row.get("time_average_stored_energy_J"),
            f"rows[{index}].time_average_stored_energy_J",
        )
        _, voltage_drop = _complex(row.get("voltage_drop_V"), f"rows[{index}].voltage_drop_V")
        p = p_complex.real
        q = q_complex.real
        w = w_complex.real
        if min(p, q, w, voltage_drop) <= 0.0:
            raise ValueError(f"rows[{index}] powers, energy, and voltage drop must be positive")

        observable_imag_relative = max(
            abs(p_complex.imag) / p,
            abs(q_complex.imag) / q,
            abs(w_complex.imag) / w,
        )
        row_metrics = {
            "index": index,
            "element_count": element_count,
            "loss_ratio_relative_error": _relative(p / q, expected_loss_ratio),
            "reactive_energy_relative_error": _relative(q, 2.0 * omega * w),
            "complex_real_relative_error": _relative(s_complex.real, p),
            "complex_imag_relative_error": _relative(s_complex.imag, q),
            "complex_magnitude_relative_error": _relative(s_abs, math.hypot(p, q)),
            "observable_imag_relative": observable_imag_relative,
            "voltage_drop_V": voltage_drop,
        }
        metrics.append(row_metrics)
        mesh_sizes.append(mesh_size)
        element_counts.append(element_count)
        powers.append(p)
        reactive_powers.append(q)
        voltage_drops.append(voltage_drop)

    max_metrics = {
        key: max(float(row[key]) for row in metrics)
        for key in (
            "loss_ratio_relative_error",
            "reactive_energy_relative_error",
            "complex_real_relative_error",
            "complex_imag_relative_error",
            "complex_magnitude_relative_error",
            "observable_imag_relative",
        )
    }
    p_convergence = _relative(powers[-2], powers[-1])
    q_convergence = _relative(reactive_powers[-2], reactive_powers[-1])
    voltage_span = (max(voltage_drops) - min(voltage_drops)) / max(voltage_drops)
    checks = {
        "three_or_more_refinements_present": len(rows) >= 3,
        "mesh_size_decreases_and_element_count_increases": all(
            right < left for left, right in zip(mesh_sizes, mesh_sizes[1:])
        )
        and all(right > left for left, right in zip(element_counts, element_counts[1:])),
        "loss_ratio_matches_sigma_over_omega_epsilon": max_metrics[
            "loss_ratio_relative_error"
        ]
        <= limits["loss_ratio_relative"],
        "reactive_power_matches_two_omega_average_energy": max_metrics[
            "reactive_energy_relative_error"
        ]
        <= limits["reactive_energy_relative"],
        "complex_power_real_part_matches_real_power": max_metrics[
            "complex_real_relative_error"
        ]
        <= limits["complex_real_relative"],
        "complex_power_imaginary_part_matches_reactive_power": max_metrics[
            "complex_imag_relative_error"
        ]
        <= limits["complex_imag_relative"],
        "complex_power_magnitude_matches_pythagorean_power": max_metrics[
            "complex_magnitude_relative_error"
        ]
        <= limits["complex_magnitude_relative"],
        "scalar_observables_are_effectively_real": max_metrics[
            "observable_imag_relative"
        ]
        <= limits["observable_imag_relative"],
        "last_two_refinements_converge": max(p_convergence, q_convergence)
        <= limits["last_pair_convergence_relative"],
        "boundary_voltage_drop_is_refinement_invariant": voltage_span
        <= limits["voltage_drop_span_relative"],
    }
    return {
        "policy": "lossy_dielectric_complex_power_refinement_gate_v1",
        "status": "ok" if all(checks.values()) else "needs_attention",
        "checks": checks,
        "issues": [name for name, ok in checks.items() if not ok],
        "limits": limits,
        "metrics": {
            "expected_sigma_over_omega_epsilon": expected_loss_ratio,
            "per_refinement": metrics,
            "maximum_errors": max_metrics,
            "last_pair_real_power_relative_change": p_convergence,
            "last_pair_reactive_power_relative_change": q_convergence,
            "voltage_drop_relative_span": voltage_span,
        },
        "lesson": (
            "A lossy-dielectric frequency-domain result must keep complex power intact: "
            "the real and imaginary parts close P and Q, its magnitude closes hypot(P,Q), "
            "and the constitutive and energy identities must survive mesh refinement."
        ),
    }
