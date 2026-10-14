"""Grounded-sphere capacitance convergence and open-boundary energy gate."""

from __future__ import annotations

import math
from collections.abc import Mapping, Sequence


_EPS0 = 8.8541878128e-12
_CASE_ORDER = ("coarse", "medium", "fine", "fine_repeat", "fine_negative")


def _finite(value: object, *, name: str) -> float:
    try:
        parsed = float(value)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"{name} must be numeric") from exc
    if not math.isfinite(parsed):
        raise ValueError(f"{name} must be finite")
    return parsed


def _mapping(value: object, *, name: str) -> Mapping[str, object]:
    if not isinstance(value, Mapping):
        raise ValueError(f"{name} must be an object")
    return value


def _rows(value: object) -> list[Mapping[str, object]]:
    if not isinstance(value, Sequence) or isinstance(value, (str, bytes)):
        raise ValueError("cases must be an array")
    rows = list(value)
    if len(rows) != len(_CASE_ORDER) or not all(
        isinstance(row, Mapping) for row in rows
    ):
        raise ValueError("exactly five case objects are required")
    return rows  # type: ignore[return-value]


def _relative_error(left: float, right: float) -> float:
    return abs(left - right) / max(abs(left), abs(right), 1.0e-300)


def _image_series_capacitance(radius_m: float, height_m: float) -> dict[str, object]:
    if radius_m <= 0.0 or height_m <= radius_m:
        raise ValueError("sphere_center_height_m must exceed a positive sphere_radius_m")
    alpha = math.acosh(height_m / radius_m)
    series_sum = 0.0
    term_count = 0
    for index in range(1, 10000):
        term = 1.0 / math.sinh(index * alpha)
        series_sum += term
        term_count = index
        if term < 1.0e-16:
            break
    capacitance = (
        4.0
        * math.pi
        * _EPS0
        * radius_m
        * math.sinh(alpha)
        * series_sum
    )
    return {
        "alpha": alpha,
        "term_count": term_count,
        "series_sum": series_sum,
        "capacitance_F": capacitance,
        "isolated_sphere_capacitance_F": 4.0 * math.pi * _EPS0 * radius_m,
    }


def _derive_case(row: Mapping[str, object], index: int) -> dict[str, object]:
    conductor = row.get("conductor")
    if not isinstance(conductor, Sequence) or isinstance(conductor, (str, bytes)):
        raise ValueError(f"cases[{index}].conductor must be [voltage, charge]")
    values = list(conductor)
    if len(values) != 2:
        raise ValueError(f"cases[{index}].conductor must contain two values")
    voltage = _finite(values[0], name=f"cases[{index}].conductor[0]")
    charge = _finite(values[1], name=f"cases[{index}].conductor[1]")
    volume_energy = _finite(
        row.get("stored_energy_J"), name=f"cases[{index}].stored_energy_J"
    )
    boundary_energy = _finite(
        row.get("mixed_boundary_energy_J"),
        name=f"cases[{index}].mixed_boundary_energy_J",
    )
    if abs(voltage) <= 0.0:
        raise ValueError(f"cases[{index}] must have nonzero voltage")
    corrected_energy = volume_energy + boundary_energy
    capacitance_charge = charge / voltage
    capacitance_volume_energy = 2.0 * volume_energy / (voltage * voltage)
    capacitance_corrected_energy = 2.0 * corrected_energy / (voltage * voltage)
    return {
        "case": str(row.get("case", "")),
        "voltage_V": voltage,
        "charge_C": charge,
        "volume_energy_J": volume_energy,
        "mixed_boundary_energy_J": boundary_energy,
        "corrected_energy_J": corrected_energy,
        "capacitance_from_charge_F": capacitance_charge,
        "capacitance_from_volume_energy_F": capacitance_volume_energy,
        "capacitance_from_corrected_energy_F": capacitance_corrected_energy,
        "charge_volume_energy_relative_error": _relative_error(
            capacitance_charge, capacitance_volume_energy
        ),
        "charge_corrected_energy_relative_error": _relative_error(
            capacitance_charge, capacitance_corrected_energy
        ),
        "node_count": int(_finite(row.get("node_count"), name=f"cases[{index}].node_count")),
        "element_count": int(
            _finite(row.get("element_count"), name=f"cases[{index}].element_count")
        ),
    }


def grounded_sphere_capacitance_convergence_gate(
    summary: Mapping[str, object],
    *,
    max_final_analytic_relative_error: float = 5.0e-4,
    max_corrected_energy_relative_error: float = 1.0e-5,
    max_replay_relative_error: float = 1.0e-13,
    max_sign_covariance_relative_error: float = 1.0e-13,
) -> dict[str, object]:
    """Validate a sphere over a ground plane with an asymptotic mixed boundary.

    The volume-field energy omits the variational contribution from
    ``eps*dV/dn + c0*V = 0``.  The total energy is therefore
    ``W_volume + 0.5*integral(c0*V^2)dS``.  The gate recomputes that identity,
    the image-sphere capacitance series, mesh convergence, deterministic replay,
    and odd-charge/even-energy covariance under voltage reversal.
    """

    if not isinstance(summary, Mapping):
        raise ValueError("summary must be an object")
    tolerances = {
        "max_final_analytic_relative_error": max_final_analytic_relative_error,
        "max_corrected_energy_relative_error": max_corrected_energy_relative_error,
        "max_replay_relative_error": max_replay_relative_error,
        "max_sign_covariance_relative_error": max_sign_covariance_relative_error,
    }
    for name, value in tolerances.items():
        parsed = _finite(value, name=name)
        if parsed < 0.0:
            raise ValueError(f"{name} must be nonnegative")
    geometry = _mapping(summary.get("problem_contract"), name="problem_contract")
    radius = _finite(geometry.get("sphere_radius_m"), name="sphere_radius_m")
    height = _finite(
        geometry.get("sphere_center_height_m"), name="sphere_center_height_m"
    )
    outer_radius = _finite(geometry.get("outer_radius_m"), name="outer_radius_m")
    c0 = _finite(
        geometry.get("open_boundary_c0_F_per_m2"),
        name="open_boundary_c0_F_per_m2",
    )
    if outer_radius <= height + radius:
        raise ValueError("outer_radius_m must enclose the sphere")
    analytic = _image_series_capacitance(radius, height)
    expected_c0 = 2.0 * _EPS0 / outer_radius
    rows = _rows(summary.get("cases"))
    derived = [_derive_case(row, index) for index, row in enumerate(rows)]

    case_names = tuple(str(row["case"]) for row in derived)
    positive = derived[:3]
    fine = derived[2]
    repeat = derived[3]
    negative = derived[4]
    analytic_errors = [
        _relative_error(
            float(row["capacitance_from_charge_F"]),
            float(analytic["capacitance_F"]),
        )
        for row in positive
    ]
    replay_errors = {
        key: _relative_error(float(fine[key]), float(repeat[key]))
        for key in (
            "capacitance_from_charge_F",
            "capacitance_from_volume_energy_F",
            "capacitance_from_corrected_energy_F",
            "volume_energy_J",
            "mixed_boundary_energy_J",
        )
    }
    sign_errors = {
        "voltage_odd": _relative_error(
            float(fine["voltage_V"]), -float(negative["voltage_V"])
        ),
        "charge_odd": _relative_error(
            float(fine["charge_C"]), -float(negative["charge_C"])
        ),
        "volume_energy_even": _relative_error(
            float(fine["volume_energy_J"]), float(negative["volume_energy_J"])
        ),
        "boundary_energy_even": _relative_error(
            float(fine["mixed_boundary_energy_J"]),
            float(negative["mixed_boundary_energy_J"]),
        ),
        "capacitance_even": _relative_error(
            float(fine["capacitance_from_charge_F"]),
            float(negative["capacitance_from_charge_F"]),
        ),
    }
    corrected_errors = [
        float(row["charge_corrected_energy_relative_error"]) for row in derived
    ]
    uncorrected_errors = [
        float(row["charge_volume_energy_relative_error"]) for row in derived
    ]
    element_counts = [int(row["element_count"]) for row in positive]
    node_counts = [int(row["node_count"]) for row in positive]
    reported_analytic = summary.get("analytic")
    if isinstance(reported_analytic, Mapping) and "capacitance_F" in reported_analytic:
        reported_analytic_error = _relative_error(
            _finite(reported_analytic["capacitance_F"], name="analytic.capacitance_F"),
            float(analytic["capacitance_F"]),
        )
    else:
        reported_analytic_error = math.inf

    checks = {
        "axisymmetric_grounded_sphere_geometry_recorded": geometry.get(
            "problem_type"
        )
        == "axisymmetric"
        and geometry.get("analysis") == "electrostatics"
        and geometry.get("length_units") == "meters"
        and _finite(
            geometry.get("ground_plane_voltage_V"), name="ground_plane_voltage_V"
        )
        == 0.0,
        "five_case_refinement_replay_sign_protocol_recorded": case_names == _CASE_ORDER,
        "asymptotic_order_two_coefficient_matches": int(
            _finite(
                geometry.get("open_boundary_asymptotic_order"),
                name="open_boundary_asymptotic_order",
            )
        )
        == 2
        and _relative_error(c0, expected_c0) <= 1.0e-8,
        "reported_image_series_matches_recomputation": reported_analytic_error
        <= 1.0e-12,
        "ground_plane_capacitance_exceeds_isolated_sphere": float(
            analytic["capacitance_F"]
        )
        > float(analytic["isolated_sphere_capacitance_F"]),
        "electrostatic_capacitances_and_energies_are_positive": all(
            int(row["node_count"]) > 0
            and int(row["element_count"]) > 0
            and float(row["capacitance_from_charge_F"]) > 0.0
            and float(row["capacitance_from_volume_energy_F"]) > 0.0
            and float(row["capacitance_from_corrected_energy_F"]) > 0.0
            and float(row["volume_energy_J"]) > 0.0
            and float(row["mixed_boundary_energy_J"]) > 0.0
            for row in derived
        ),
        "realized_mesh_refines_monotonically": node_counts[0]
        < node_counts[1]
        < node_counts[2]
        and element_counts[0] < element_counts[1] < element_counts[2],
        "charge_capacitance_converges_to_image_series": analytic_errors[0]
        > analytic_errors[1]
        > analytic_errors[2]
        and analytic_errors[2] <= float(max_final_analytic_relative_error),
        "volume_energy_alone_reveals_boundary_omission": min(uncorrected_errors)
        >= 1.0e-3,
        "mixed_boundary_energy_restores_charge_energy_identity": max(
            corrected_errors
        )
        <= float(max_corrected_energy_relative_error)
        and max(corrected_errors) <= 0.01 * min(uncorrected_errors),
        "fine_replay_is_deterministic": max(replay_errors.values())
        <= float(max_replay_relative_error)
        and int(fine["node_count"]) == int(repeat["node_count"])
        and int(fine["element_count"]) == int(repeat["element_count"]),
        "voltage_reversal_has_odd_charge_and_even_energy": max(sign_errors.values())
        <= float(max_sign_covariance_relative_error)
        and int(fine["node_count"]) == int(negative["node_count"])
        and int(fine["element_count"]) == int(negative["element_count"]),
    }
    issues = [name for name, passed in checks.items() if not passed]
    return {
        "policy": "grounded_sphere_capacitance_convergence_gate_v1",
        "status": "ok" if not issues else "needs_attention",
        "checks": checks,
        "issues": issues,
        "analytic": analytic,
        "derived_cases": derived,
        "metrics": {
            "analytic_relative_errors": analytic_errors,
            "maximum_corrected_energy_relative_error": max(corrected_errors),
            "minimum_uncorrected_energy_relative_error": min(uncorrected_errors),
            "maximum_replay_relative_error": max(replay_errors.values()),
            "maximum_sign_covariance_relative_error": max(sign_errors.values()),
            "open_boundary_coefficient_relative_error": _relative_error(c0, expected_c0),
            "reported_analytic_relative_error": reported_analytic_error,
        },
        "lesson": (
            "For an axisymmetric sphere above a grounded plane, compare terminal charge "
            "with the image-sphere series. A mixed outer boundary contributes "
            "0.5*integral(c0*V^2)dS to the variational energy; omitting it can make an "
            "otherwise accurate charge result appear inconsistent with field energy."
        ),
    }
