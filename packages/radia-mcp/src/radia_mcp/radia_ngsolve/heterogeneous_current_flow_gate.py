"""Solver-neutral gate for heterogeneous current-flow P1 reintegration."""

from __future__ import annotations

import math
from collections.abc import Mapping, Sequence
from typing import Any


_CASE_ORDER = ("coarse", "medium", "fine", "fine_repeat", "fine_negative")
_ODD_POINT_FIELDS = frozenset(
    {"V", "Jx", "Jy", "Ex", "Ey", "Jdx", "Jdy", "Jcx", "Jcy"}
)
_DEFAULT_LIMITS = {
    "identity_relative": 1.0e-8,
    "reintegration_relative": 1.0e-8,
    "material_partition_relative": 1.0e-10,
    "last_pair_convergence_relative": 2.0e-3,
    "exact_replay_relative": 1.0e-12,
    "sign_covariance_relative": 1.0e-10,
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


def _integer(value: object, name: str, *, positive: bool = False) -> int:
    parsed = _finite(value, name, positive=positive)
    integer = int(parsed)
    if parsed != integer:
        raise ValueError(f"{name} must be an integer")
    return integer


def _mapping(value: object, name: str) -> Mapping[str, object]:
    if not isinstance(value, Mapping):
        raise ValueError(f"{name} must be an object")
    return value


def _sequence(value: object, name: str) -> list[object]:
    if not isinstance(value, Sequence) or isinstance(value, (str, bytes)):
        raise ValueError(f"{name} must be an array")
    return list(value)


def _relative(actual: float | complex, reference: float | complex) -> float:
    return abs(actual - reference) / max(abs(reference), abs(actual), 1.0e-300)


def _complex(value: object, name: str) -> complex:
    if not isinstance(value, Mapping):
        return complex(_finite(value, name), 0.0)
    encoded = _mapping(value, name)
    number = complex(
        _finite(encoded.get("real"), f"{name}.real"),
        _finite(encoded.get("imag"), f"{name}.imag"),
    )
    magnitude = _finite(encoded.get("abs"), f"{name}.abs")
    if magnitude < 0.0:
        raise ValueError(f"{name}.abs must be nonnegative")
    if _relative(abs(number), magnitude) > 1.0e-12:
        raise ValueError(f"{name}.abs is inconsistent with real and imag")
    return number


def _case_metrics(
    row: Mapping[str, object],
    index: int,
    *,
    omega: float,
    depth: float,
) -> dict[str, Any]:
    prefix = f"cases[{index}]"
    mesh_size = _finite(
        row.get("mesh_size_mm"), f"{prefix}.mesh_size_mm", positive=True
    )
    element_count = _integer(
        row.get("element_count"), f"{prefix}.element_count", positive=True
    )
    real_power = _complex(row.get("real_power_W"), f"{prefix}.real_power_W")
    apparent_power = _complex(
        row.get("apparent_power_VA"), f"{prefix}.apparent_power_VA"
    )
    energy = _complex(
        row.get("time_average_stored_energy_J"),
        f"{prefix}.time_average_stored_energy_J",
    )
    area = _complex(row.get("area_m2"), f"{prefix}.area_m2")
    volume = _complex(row.get("volume_m3"), f"{prefix}.volume_m3")
    hv_voltage = _complex(row.get("hv_voltage_V"), f"{prefix}.hv_voltage_V")
    high_voltage = _finite(row.get("high_voltage_V"), f"{prefix}.high_voltage_V")
    if min(real_power.real, apparent_power.imag, energy.real, area.real, volume.real) <= 0.0:
        raise ValueError(f"{prefix} power, energy, area, and volume must be positive")

    local = _mapping(
        row.get("maximum_local_identity_errors"),
        f"{prefix}.maximum_local_identity_errors",
    )
    local_errors = {
        key: _finite(local.get(key), f"{prefix}.local_errors.{key}")
        for key in (
            "total_current_split_relative_error",
            "complex_conductivity_relative_error",
            "conduction_current_relative_error",
            "displacement_current_relative_error",
        )
    }
    if min(local_errors.values()) < 0.0:
        raise ValueError(f"{prefix} local identity errors must be nonnegative")

    independent = _mapping(
        row.get("independent_anc_reintegration"),
        f"{prefix}.independent_anc_reintegration",
    )
    independent_elements = _integer(
        independent.get("element_count"), f"{prefix}.anc.element_count", positive=True
    )
    independent_area = _finite(
        independent.get("total_area_m2"), f"{prefix}.anc.total_area_m2", positive=True
    )
    independent_power = _complex(
        independent.get("total_complex_power_VA"),
        f"{prefix}.anc.total_complex_power_VA",
    )
    independent_energy = _finite(
        independent.get("total_energy_J"),
        f"{prefix}.anc.total_energy_J",
        positive=True,
    )
    independent_two_omega_energy = _finite(
        independent.get("two_omega_energy_var"),
        f"{prefix}.anc.two_omega_energy_var",
        positive=True,
    )
    material_rows = _mapping(
        independent.get("material_rows"), f"{prefix}.anc.material_rows"
    )
    if len(material_rows) < 2:
        raise ValueError(f"{prefix} must contain at least two material partitions")

    material_counts: list[int] = []
    material_areas: list[float] = []
    material_powers: list[complex] = []
    material_energies: list[float] = []
    for material_name, material_value in material_rows.items():
        material = _mapping(material_value, f"{prefix}.materials.{material_name}")
        material_counts.append(
            _integer(
                material.get("element_count"),
                f"{prefix}.materials.{material_name}.element_count",
                positive=True,
            )
        )
        material_areas.append(
            _finite(
                material.get("area_m2"),
                f"{prefix}.materials.{material_name}.area_m2",
                positive=True,
            )
        )
        material_powers.append(
            _complex(
                material.get("complex_power_VA"),
                f"{prefix}.materials.{material_name}.complex_power_VA",
            )
        )
        material_energies.append(
            _finite(
                material.get("energy_J"),
                f"{prefix}.materials.{material_name}.energy_J",
            )
        )

    postprocessor_errors = {
        "real_power_is_complex_power_real": _relative(real_power.real, apparent_power.real),
        "reactive_power_is_two_omega_energy": _relative(
            apparent_power.imag, 2.0 * omega * energy.real
        ),
        "area_times_depth_is_volume": _relative(area.real * depth, volume.real),
        "scalar_observables_are_real": max(
            _relative(real_power.imag, 0.0),
            _relative(energy.imag, 0.0),
            _relative(area.imag, 0.0),
            _relative(volume.imag, 0.0),
        ),
    }
    reintegration_errors = {
        "element_count": 0.0 if independent_elements == element_count else math.inf,
        "real_power": _relative(independent_power.real, real_power.real),
        "reactive_power": _relative(independent_power.imag, apparent_power.imag),
        "stored_energy": _relative(independent_energy, energy.real),
        "area": _relative(independent_area, area.real),
        "energy_identity": _relative(
            independent_power.imag, 2.0 * omega * independent_energy
        ),
        "stored_two_omega_energy": _relative(
            independent_two_omega_energy, 2.0 * omega * independent_energy
        ),
    }
    material_power = sum(material_powers, start=0j)
    partition_errors = {
        "element_count": 0.0 if sum(material_counts) == element_count else math.inf,
        "area": _relative(sum(material_areas), independent_area),
        "complex_power": _relative(material_power, independent_power),
        "energy": _relative(sum(material_energies), independent_energy),
    }
    return {
        "name": str(row.get("case", "")),
        "mesh_size_mm": mesh_size,
        "element_count": element_count,
        "high_voltage_V": high_voltage,
        "real_power": real_power,
        "apparent_power": apparent_power,
        "energy": energy,
        "hv_voltage": hv_voltage,
        "maximum_local_error": max(local_errors.values()),
        "maximum_postprocessor_error": max(postprocessor_errors.values()),
        "maximum_reintegration_error": max(reintegration_errors.values()),
        "maximum_partition_error": max(partition_errors.values()),
        "material_count": len(material_rows),
        "point_rows": row.get("point_rows"),
        "details": {
            "local_identity_errors": local_errors,
            "postprocessor_errors": postprocessor_errors,
            "reintegration_errors": reintegration_errors,
            "material_partition_errors": partition_errors,
        },
    }


def _replay_error(left: Mapping[str, Any], right: Mapping[str, Any]) -> float:
    return max(
        _relative(left["real_power"], right["real_power"]),
        _relative(left["apparent_power"], right["apparent_power"]),
        _relative(left["energy"], right["energy"]),
        _relative(left["hv_voltage"], right["hv_voltage"]),
    )


def _point_sign_error(
    positive_rows: object,
    negative_rows: object,
    point_value_order: Sequence[str],
) -> float:
    positive_sets = (
        [positive_rows]
        if isinstance(positive_rows, Mapping)
        else _sequence(positive_rows, "fine.point_rows")
    )
    negative_sets = (
        [negative_rows]
        if isinstance(negative_rows, Mapping)
        else _sequence(negative_rows, "fine_negative.point_rows")
    )
    if len(positive_sets) != len(negative_sets) or not positive_sets:
        raise ValueError("fine and fine_negative point_rows must have equal nonzero length")
    maximum = 0.0
    for set_index, (positive_set, negative_set) in enumerate(
        zip(positive_sets, negative_sets)
    ):
        positive_map = _mapping(positive_set, f"fine.point_rows[{set_index}]")
        negative_map = _mapping(negative_set, f"fine_negative.point_rows[{set_index}]")
        if set(positive_map) != set(negative_map):
            raise ValueError("fine and fine_negative point material names must match")
        for material_name, positive_value in positive_map.items():
            positive = _mapping(positive_value, f"fine.{material_name}")
            negative = _mapping(negative_map[material_name], f"fine_negative.{material_name}")
            if positive.get("point_mm") != negative.get("point_mm"):
                raise ValueError("fine and fine_negative point coordinates must match")
            positive_values = _sequence(positive.get("values"), f"fine.{material_name}.values")
            negative_values = _sequence(negative.get("values"), f"fine_negative.{material_name}.values")
            if len(positive_values) != len(point_value_order) or len(negative_values) != len(point_value_order):
                raise ValueError("point value count must match point_value_order")
            for value_index, field_name in enumerate(point_value_order):
                if field_name not in _ODD_POINT_FIELDS:
                    continue
                left = _complex(
                    positive_values[value_index], f"fine.{material_name}.{field_name}"
                )
                right = _complex(
                    negative_values[value_index],
                    f"fine_negative.{material_name}.{field_name}",
                )
                maximum = max(maximum, _relative(left, -right))
    return maximum


def heterogeneous_current_flow_p1_reintegration_gate(
    summary: Mapping[str, object],
) -> dict[str, Any]:
    """Gate heterogeneous current flow using independent P1 triangle evidence."""

    if not isinstance(summary, Mapping):
        raise ValueError("summary must be an object")
    frequency = _finite(summary.get("frequency_Hz"), "frequency_Hz", positive=True)
    depth = _finite(summary.get("depth_m"), "depth_m", positive=True)
    omega = 2.0 * math.pi * frequency
    overrides = summary.get("gate_tolerances", {})
    if not isinstance(overrides, Mapping):
        raise ValueError("gate_tolerances must be an object")
    limits = {
        key: _finite(overrides.get(key, default), key, positive=True)
        for key, default in _DEFAULT_LIMITS.items()
    }

    postprocess = _mapping(summary.get("postprocess_contract"), "postprocess_contract")
    point_value_order = [
        str(value)
        for value in _sequence(
            postprocess.get("point_value_order"), "point_value_order"
        )
    ]
    if not _ODD_POINT_FIELDS.issubset(point_value_order):
        raise ValueError("point_value_order is missing odd field observables")
    anc = _mapping(summary.get("anc_contract"), "anc_contract")
    p1_contract = (
        anc.get("element_order") == "P1_triangle"
        and isinstance(anc.get("material_resolution"), str)
        and bool(anc.get("material_resolution"))
        and isinstance(anc.get("independent_power_identity"), str)
        and bool(anc.get("independent_power_identity"))
    )

    rows = _sequence(summary.get("cases"), "cases")
    if len(rows) != len(_CASE_ORDER) or any(not isinstance(row, Mapping) for row in rows):
        raise ValueError("cases must contain exactly five result objects")
    case_metrics = [
        _case_metrics(row, index, omega=omega, depth=depth)
        for index, row in enumerate(rows)
    ]
    if tuple(metric["name"] for metric in case_metrics) != _CASE_ORDER:
        raise ValueError(f"cases must be ordered as {', '.join(_CASE_ORDER)}")

    coarse, medium, fine, repeat, negative = case_metrics
    convergence = max(
        _relative(medium["real_power"], fine["real_power"]),
        _relative(medium["apparent_power"].imag, fine["apparent_power"].imag),
    )
    replay_error = _replay_error(fine, repeat)
    sign_error = max(
        _relative(fine["high_voltage_V"], -negative["high_voltage_V"]),
        _relative(fine["hv_voltage"], -negative["hv_voltage"]),
        _relative(fine["real_power"], negative["real_power"]),
        _relative(fine["apparent_power"], negative["apparent_power"]),
        _relative(fine["energy"], negative["energy"]),
        _point_sign_error(
            fine["point_rows"], negative["point_rows"], point_value_order
        ),
    )
    maximum_local = max(metric["maximum_local_error"] for metric in case_metrics)
    maximum_postprocessor = max(
        metric["maximum_postprocessor_error"] for metric in case_metrics
    )
    maximum_reintegration = max(
        metric["maximum_reintegration_error"] for metric in case_metrics
    )
    maximum_partition = max(
        metric["maximum_partition_error"] for metric in case_metrics
    )
    checks = {
        "p1_triangle_reintegration_contract_is_explicit": p1_contract,
        "three_refinement_levels_are_monotone": all(
            right["mesh_size_mm"] < left["mesh_size_mm"]
            and right["element_count"] > left["element_count"]
            for left, right in zip((coarse, medium), (medium, fine))
        ),
        "heterogeneous_material_partitions_are_present": all(
            metric["material_count"] >= 2 for metric in case_metrics
        ),
        "pointwise_current_identities_close": maximum_local
        <= limits["identity_relative"],
        "global_power_energy_and_geometry_identities_close": maximum_postprocessor
        <= limits["identity_relative"],
        "independent_p1_reintegration_closes": maximum_reintegration
        <= limits["reintegration_relative"],
        "material_partition_sums_close": maximum_partition
        <= limits["material_partition_relative"],
        "medium_to_fine_power_converges": convergence
        <= limits["last_pair_convergence_relative"],
        "fine_replay_is_exact": replay_error <= limits["exact_replay_relative"],
        "voltage_and_fields_are_odd_while_quadratic_outputs_are_even": sign_error
        <= limits["sign_covariance_relative"],
    }
    return {
        "policy": "heterogeneous_current_flow_p1_reintegration_gate_v1",
        "status": "ok" if all(checks.values()) else "needs_attention",
        "checks": checks,
        "issues": [name for name, ok in checks.items() if not ok],
        "limits": limits,
        "metrics": {
            "element_counts": [metric["element_count"] for metric in case_metrics],
            "maximum_local_identity_relative_error": maximum_local,
            "maximum_postprocessor_identity_relative_error": maximum_postprocessor,
            "maximum_independent_p1_reintegration_relative_error": maximum_reintegration,
            "maximum_material_partition_relative_error": maximum_partition,
            "medium_to_fine_power_relative_change": convergence,
            "fine_replay_relative_error": replay_error,
            "sign_covariance_relative_error": sign_error,
            "per_case": [
                {
                    "case": metric["name"],
                    "mesh_size_mm": metric["mesh_size_mm"],
                    "element_count": metric["element_count"],
                    "material_count": metric["material_count"],
                    **metric["details"],
                }
                for metric in case_metrics
            ],
        },
        "lesson": (
            "A heterogeneous current-flow result needs pointwise J=Jc+Jd, global "
            "P+jQ and energy closure, independent P1 material reintegration, mesh "
            "convergence, deterministic replay, and odd-field/even-power covariance."
        ),
    }
