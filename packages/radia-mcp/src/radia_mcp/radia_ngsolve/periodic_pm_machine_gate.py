"""Solver-neutral periodicity and replay gate for unwrapped PM machines."""

from __future__ import annotations

import math
from collections.abc import Mapping


def _relative(left: float, right: float) -> float:
    return abs(left - right) / max(abs(left), abs(right), 1.0e-300)


def periodic_unwrapped_pm_machine_replay_gate(
    summary: Mapping[str, object],
    *,
    maximum_field_symmetry_relative_error: float = 0.05,
    maximum_energy_replay_relative_error: float = 1.0e-3,
    maximum_field_replay_relative_error: float = 5.0e-3,
    maximum_mesh_cardinality_relative_difference: float = 0.05,
) -> dict[str, object]:
    """Gate topology-aware half-turn symmetry and observable replay stability."""
    if not isinstance(summary, Mapping):
        raise TypeError("summary must be a mapping")
    tolerances = {
        "maximum_field_symmetry_relative_error": float(maximum_field_symmetry_relative_error),
        "maximum_energy_replay_relative_error": float(maximum_energy_replay_relative_error),
        "maximum_field_replay_relative_error": float(maximum_field_replay_relative_error),
        "maximum_mesh_cardinality_relative_difference": float(
            maximum_mesh_cardinality_relative_difference
        ),
    }
    if not all(math.isfinite(value) and value >= 0.0 for value in tolerances.values()):
        raise ValueError("tolerances must be finite and nonnegative")

    machine = summary.get("machine")
    boundary = summary.get("periodic_boundary")
    records = summary.get("runs")
    if not isinstance(machine, Mapping) or not isinstance(boundary, Mapping):
        raise ValueError("machine and periodic_boundary must be mappings")
    if not isinstance(records, list) or len(records) != 2:
        raise ValueError("runs must contain exactly two records")

    slot_count = int(machine.get("slot_count", 0))
    pole_count = int(machine.get("pole_count", 0))
    circumference = float(machine.get("circumference", math.nan))
    slot_pitch = float(machine.get("slot_pitch", math.nan))
    pole_pitch = float(machine.get("pole_pitch", math.nan))
    if not all(math.isfinite(value) and value > 0.0 for value in (circumference, slot_pitch, pole_pitch)):
        raise ValueError("circumference and pitches must be finite and positive")

    parsed = []
    for index, record in enumerate(records):
        if not isinstance(record, Mapping):
            raise ValueError(f"runs[{index}] must be a mapping")
        values = {
            "energy": float(record.get("energy", math.nan)),
            "coenergy": float(record.get("coenergy", math.nan)),
            "normal_flux_rms": float(record.get("normal_flux_rms", math.nan)),
            "normal_flux_mean_relative": float(
                record.get("normal_flux_mean_relative", math.nan)
            ),
            "half_turn_antiperiodicity_relative_error": float(
                record.get("half_turn_antiperiodicity_relative_error", math.nan)
            ),
            "finite_profile_coverage": float(record.get("finite_profile_coverage", math.nan)),
            "antiperiodic_pair_coverage": float(
                record.get("antiperiodic_pair_coverage", math.nan)
            ),
            "node_count": float(record.get("node_count", math.nan)),
            "element_count": float(record.get("element_count", math.nan)),
        }
        if not all(math.isfinite(value) for value in values.values()):
            raise ValueError(f"runs[{index}] contains non-finite observables")
        parsed.append(values)

    first, second = parsed
    replay_errors = {
        "energy": _relative(first["energy"], second["energy"]),
        "coenergy": _relative(first["coenergy"], second["coenergy"]),
        "normal_flux_rms": _relative(
            first["normal_flux_rms"], second["normal_flux_rms"]
        ),
        "normal_flux_mean_relative": _relative(
            first["normal_flux_mean_relative"], second["normal_flux_mean_relative"]
        ),
        "node_count": _relative(first["node_count"], second["node_count"]),
        "element_count": _relative(first["element_count"], second["element_count"]),
    }
    checks = {
        "pitch_closure": _relative(slot_count * slot_pitch, circumference) <= 1.0e-12
        and _relative(pole_count * pole_pitch, circumference) <= 1.0e-12,
        "half_turn_slot_geometry_repeats": slot_count > 0 and slot_count % 2 == 0,
        "half_turn_magnetization_reverses": pole_count > 0
        and pole_count % 2 == 0
        and (pole_count // 2) % 2 == 1,
        "symmetry_shift_is_topology_aware": machine.get("field_symmetry_shift")
        == "half_circumference",
        "periodic_boundary_pair_nonempty": int(boundary.get("master_count", 0)) > 0
        and int(boundary.get("master_count", 0))
        == int(boundary.get("slave_count", -1)),
        "p1_triangle_meshes_recorded": all(
            str(record.get("element_type") or "") == "TL3"
            and float(record.get("node_count", 0.0)) > 100.0
            and float(record.get("element_count", 0.0)) > 100.0
            for record in records
        ),
        "positive_energy_and_coenergy": all(
            record["energy"] > 0.0 and record["coenergy"] > 0.0 for record in parsed
        ),
        "airgap_profiles_well_covered": all(
            record["finite_profile_coverage"] >= 0.95
            and record["antiperiodic_pair_coverage"] >= 0.90
            for record in parsed
        ),
        "low_net_normal_flux": all(
            record["normal_flux_mean_relative"]
            <= tolerances["maximum_field_symmetry_relative_error"]
            for record in parsed
        ),
        "half_turn_antiperiodicity_closes": all(
            record["half_turn_antiperiodicity_relative_error"]
            <= tolerances["maximum_field_symmetry_relative_error"]
            for record in parsed
        ),
        "energy_observables_replay": max(replay_errors["energy"], replay_errors["coenergy"])
        <= tolerances["maximum_energy_replay_relative_error"],
        "field_observables_replay": max(
            replay_errors["normal_flux_rms"],
            replay_errors["normal_flux_mean_relative"],
        )
        <= tolerances["maximum_field_replay_relative_error"],
        "mesh_cardinality_remains_same_scale": max(
            replay_errors["node_count"], replay_errors["element_count"]
        )
        <= tolerances["maximum_mesh_cardinality_relative_difference"],
    }
    issues = [name for name, ok in checks.items() if not ok]
    return {
        "policy": "periodic_unwrapped_pm_machine_replay_gate_v1",
        "status": "ok" if not issues else "needs_attention",
        "checks": checks,
        "issues": issues,
        "metrics": {"replay_relative_errors": replay_errors},
        "tolerances": tolerances,
        "lesson": (
            "Field symmetry must combine stator-slot repetition and rotor-pole polarity. "
            "A one-pole shift is not generally a symmetry of the complete machine. When "
            "meshing is nondeterministic, gate reproducibility on physical observables and "
            "treat mesh cardinality only as a same-scale diagnostic."
        ),
    }
