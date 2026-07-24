"""Strict artifact gate for the public TEAM Workshop Problem 36 benchmark."""

from __future__ import annotations

import hashlib
import json
import math
from collections.abc import Mapping, Sequence
from typing import Any


POLICY = "radia_ih_team36_axisymmetric_gate_v1"
SCHEMA = "radia_ih_team36_axisymmetric_v1"
SOURCE_URL = "https://www.compumag.org/wp/wp-content/uploads/2021/07/problem-36.pdf"

EXPECTED_GEOMETRY_M = {
    "billet_length": 1.0,
    "billet_radius": 0.03,
    "coil_length": 1.0,
    "coil_inner_radius": 0.048,
    "turn_axial_length": 0.04,
    "turn_radial_width": 0.02,
    "turn_conductor_thickness": 0.003,
}
EXPECTED_EXCITATION = {
    "frequency_hz": 2000.0,
    "current_rms_a": 3500.0,
    "duration_s": 250.0,
}
REFERENCE_OBSERVABLES = {
    "axis_temperature_c_at_250_s": "axis_temperature_c",
    "surface_temperature_c_at_250_s": "surface_temperature_c",
    "maximum_temperature_c_at_250_s": "maximum_temperature_c",
    "induced_power_w_at_250_s": "induced_power_w",
}


def _finite(value: object) -> bool:
    try:
        return math.isfinite(float(value))
    except (TypeError, ValueError):
        return False


def _close(value: object, expected: float, *, rel: float = 1.0e-12) -> bool:
    return _finite(value) and math.isclose(
        float(value), expected, rel_tol=rel, abs_tol=max(1.0e-15, rel * abs(expected))
    )


def _sha256(value: object) -> bool:
    text = str(value).lower()
    return len(text) == 64 and all(char in "0123456789abcdef" for char in text)


def _git_digest(value: object) -> bool:
    text = str(value).lower()
    return (
        len(text) in {40, 64}
        and any(char != "0" for char in text)
        and all(char in "0123456789abcdef" for char in text)
    )


def _mapping(value: object) -> Mapping[str, object]:
    return value if isinstance(value, Mapping) else {}


def _sequence(value: object) -> Sequence[object]:
    if isinstance(value, Sequence) and not isinstance(value, (str, bytes, bytearray)):
        return value
    return ()


def _same_identity(left: Mapping[str, object], right: Mapping[str, object]) -> bool:
    fields = (
        "geometry_sha256",
        "material_tables_sha256",
        "excitation_sha256",
        "coordinate_system",
    )
    return all(left.get(field) == right.get(field) for field in fields)


def _history_checks(history: Sequence[object]) -> dict[str, bool]:
    rows = [row for row in history if isinstance(row, Mapping)]
    times = [float(row["time_s"]) for row in rows if _finite(row.get("time_s"))]
    complete_rows = len(rows) == len(history) and len(rows) >= 2
    finite_observables = complete_rows and all(
        _finite(row.get("axis_temperature_c"))
        and _finite(row.get("surface_temperature_c"))
        and _finite(row.get("maximum_temperature_c"))
        for row in rows
    )
    increasing = (
        len(times) == len(rows)
        and all(right > left for left, right in zip(times, times[1:]))
    )
    positive_power = complete_rows and all(
        _finite(row.get("induced_power_w")) and float(row["induced_power_w"]) > 0.0
        for row in rows[1:]
    )
    em_converged = complete_rows and all(bool(row.get("em_converged")) for row in rows[1:])
    thermal_converged = complete_rows and all(
        bool(row.get("thermal_converged")) for row in rows[1:]
    )
    return {
        "history_has_initial_and_transient_rows": complete_rows,
        "history_starts_at_zero": bool(times) and _close(times[0], 0.0),
        "history_ends_at_250_s": bool(times) and _close(times[-1], 250.0),
        "history_time_is_strictly_increasing": increasing,
        "history_observables_are_finite": finite_observables,
        "induced_power_is_positive_after_start": positive_power,
        "all_em_steps_converged": em_converged,
        "all_thermal_steps_converged": thermal_converged,
    }


def _reference_checks(
    artifact: Mapping[str, object], reference: Mapping[str, object]
) -> dict[str, bool]:
    if not reference:
        return {
            "cross_reference_supplied": False,
            "cross_reference_identity_matches": False,
            "cross_reference_observables_match": False,
        }
    artifact_identity = _mapping(artifact.get("identity"))
    reference_identity = _mapping(reference.get("identity"))
    comparisons = _sequence(reference.get("comparisons"))
    comparison_ok = bool(comparisons)
    for item in comparisons:
        row = _mapping(item)
        observable = str(row.get("observable", ""))
        artifact_value = _artifact_observable(artifact, observable)
        absolute_scale = row.get("absolute_scale", 1.0e-12)
        if not (
            artifact_value is not None
            and _finite(row.get("reference_value"))
            and _finite(row.get("relative_tolerance"))
            and _finite(absolute_scale)
        ):
            comparison_ok = False
            break
        expected = float(row["reference_value"])
        tolerance = float(row["relative_tolerance"])
        scale_floor = float(absolute_scale)
        reported_value = row.get("radia_value")
        reported_matches = reported_value is None or (
            _finite(reported_value)
            and _close(reported_value, artifact_value, rel=1.0e-12)
        )
        scale = max(abs(expected), scale_floor)
        relative_error = abs(artifact_value - expected) / scale
        if (
            tolerance < 0.0
            or scale_floor <= 0.0
            or not reported_matches
            or relative_error > tolerance
        ):
            comparison_ok = False
            break
    return {
        "cross_reference_supplied": True,
        "cross_reference_identity_matches": _same_identity(
            artifact_identity, reference_identity
        ),
        "cross_reference_observables_match": comparison_ok,
    }


def _artifact_observable(
    artifact: Mapping[str, object], observable: str
) -> float | None:
    field = REFERENCE_OBSERVABLES.get(observable)
    if field is None:
        return None
    rows = [
        row
        for row in _sequence(artifact.get("history"))
        if isinstance(row, Mapping) and _close(row.get("time_s"), 250.0)
    ]
    if len(rows) != 1 or not _finite(rows[0].get(field)):
        return None
    return float(rows[0][field])


def team36_contract() -> dict[str, Any]:
    """Return the immutable public benchmark contract used by the gate."""

    return {
        "policy": POLICY,
        "artifact_schema": SCHEMA,
        "benchmark_source": SOURCE_URL,
        "geometry_m": dict(EXPECTED_GEOMETRY_M),
        "excitation": dict(EXPECTED_EXCITATION),
        "required_physics": [
            "2d_axisymmetric_time_harmonic_em",
            "transient_thermal_250_s",
            "mu_of_field_and_temperature",
            "electrical_resistivity_of_temperature",
            "thermal_conductivity_of_temperature",
            "heat_capacity_of_temperature",
            "convection_and_radiation",
        ],
        "required_coupling": [
            "noncoincident_em_and_thermal_meshes",
            "temperature_to_em_mapping",
            "joule_power_to_thermal_conservative_mapping",
        ],
        "acceptance": {
            "power_mapping_relative_error_max": 0.02,
            "cross_validation_requires_independent_reference": True,
            "reference_observables": sorted(REFERENCE_OBSERVABLES),
        },
    }


def evaluate_team36_artifact(
    artifact: Mapping[str, object],
    *,
    reference: Mapping[str, object] | None = None,
) -> dict[str, Any]:
    """Evaluate one saved TEAM 36 result without reading arbitrary local paths."""

    geometry = _mapping(artifact.get("geometry_m"))
    excitation = _mapping(artifact.get("excitation"))
    material = _mapping(artifact.get("material_model"))
    meshes = _mapping(artifact.get("meshes"))
    em_mesh = _mapping(meshes.get("electromagnetic"))
    thermal_mesh = _mapping(meshes.get("thermal"))
    coupling = _mapping(artifact.get("coupling"))
    temperature_map = _mapping(coupling.get("temperature_to_em"))
    power_map = _mapping(coupling.get("joule_power_to_thermal"))
    provenance = _mapping(artifact.get("provenance"))
    timing = _mapping(artifact.get("timing_s"))

    checks: dict[str, bool] = {
        "result_starts_with_version_date_and_host": list(artifact)[:3]
        == ["radia_version", "executed_at_utc", "host"]
        and artifact.get("radia_version") == provenance.get("radia_version")
        and artifact.get("executed_at_utc") == provenance.get("executed_at_utc")
        and artifact.get("host") == provenance.get("host"),
        "schema_matches": artifact.get("artifact_schema") == SCHEMA,
        "benchmark_source_is_public_original": artifact.get("benchmark_source") == SOURCE_URL,
        "axisymmetric_coordinate_system": artifact.get("coordinate_system") == "axisymmetric_r_z",
        "geometry_matches_original": all(
            _close(geometry.get(name), value)
            for name, value in EXPECTED_GEOMETRY_M.items()
        ),
        "excitation_matches_original": all(
            _close(excitation.get(name), value)
            for name, value in EXPECTED_EXCITATION.items()
        ),
        "material_tables_are_complete": all(
            int(material.get(name, 0)) >= minimum
            for name, minimum in {
                "resistivity_point_count": 16,
                "mu20_point_count": 20,
                "conductivity_point_count": 17,
                "heat_capacity_point_count": 30,
            }.items()
        )
        and _close(material.get("curie_temperature_c"), 770.0)
        and _close(material.get("transition_width_c"), 20.0),
        "mesh_hashes_are_present": _sha256(em_mesh.get("topology_sha256"))
        and _sha256(thermal_mesh.get("topology_sha256")),
        "meshes_are_first_order_triangles": em_mesh.get("element_kinds") == ["triangle"]
        and thermal_mesh.get("element_kinds") == ["triangle"]
        and int(em_mesh.get("element_order", 0)) == 1
        and int(thermal_mesh.get("element_order", 0)) == 1
        and int(em_mesh.get("billet_skin_layer_count", 0)) >= 4,
        "meshes_are_noncoincident": em_mesh.get("topology_sha256")
        != thermal_mesh.get("topology_sha256")
        and (
            em_mesh.get("vertex_count") != thermal_mesh.get("vertex_count")
            or em_mesh.get("element_count") != thermal_mesh.get("element_count")
        ),
        "temperature_mapping_is_bidirectional_evidence": temperature_map.get("source_mesh_sha256")
        == thermal_mesh.get("topology_sha256")
        and temperature_map.get("target_mesh_sha256") == em_mesh.get("topology_sha256")
        and int(temperature_map.get("sample_count", 0)) > 0
        and int(temperature_map.get("outside_count", -1)) == 0,
        "joule_mapping_is_conservative": power_map.get("source_mesh_sha256")
        == em_mesh.get("topology_sha256")
        and power_map.get("target_mesh_sha256") == thermal_mesh.get("topology_sha256")
        and int(power_map.get("sample_count", 0)) > 0
        and _finite(power_map.get("maximum_relative_error"))
        and float(power_map["maximum_relative_error"]) <= 0.02
        and _finite(power_map.get("maximum_scale_deviation"))
        and float(power_map["maximum_scale_deviation"]) <= 0.25,
        "run_timestamp_and_versions_present": bool(provenance.get("executed_at_utc"))
        and bool(provenance.get("host"))
        and bool(provenance.get("radia_version"))
        and bool(provenance.get("ngsolve_version"))
        and _git_digest(provenance.get("git_commit")),
        "four_major_timing_components_present": all(
            _finite(timing.get(name)) and float(timing[name]) >= 0.0
            for name in ("mesh", "electromagnetic", "mapping", "thermal")
        ),
    }
    checks.update(_history_checks(_sequence(artifact.get("history"))))
    reference_checks = _reference_checks(artifact, reference or {})
    checks.update(reference_checks)

    solver_checks = {
        name: result
        for name, result in checks.items()
        if not name.startswith("cross_reference_")
    }
    accepted_for_solver_execution = all(solver_checks.values())
    accepted_for_cross_validation = accepted_for_solver_execution and all(
        reference_checks.values()
    )
    return {
        "policy": POLICY,
        "accepted_for_solver_execution": accepted_for_solver_execution,
        "accepted_for_cross_validation": accepted_for_cross_validation,
        "accepted_for_mcp_learning": accepted_for_cross_validation,
        "checks": checks,
        "failed_checks": [name for name, result in checks.items() if not result],
        "next_action": (
            "Artifact and independent reference agree; promote both as a regression pair."
            if accepted_for_cross_validation
            else "Supply an identity-matched independent reference after all solver checks pass."
            if accepted_for_solver_execution
            else "Repair the failed solver or artifact checks before comparing solvers."
        ),
    }


def evaluate_team36_json(artifact_json: str, reference_json: str = "") -> dict[str, Any]:
    """JSON-only MCP boundary; deliberately does not accept filesystem paths."""

    artifact = json.loads(artifact_json)
    reference = json.loads(reference_json) if reference_json.strip() else None
    if not isinstance(artifact, Mapping):
        raise ValueError("artifact_json must encode an object")
    if reference is not None and not isinstance(reference, Mapping):
        raise ValueError("reference_json must encode an object")
    return evaluate_team36_artifact(artifact, reference=reference)


def identity_digest(value: Mapping[str, object]) -> str:
    """Canonical SHA-256 helper shared by the validation runner."""

    payload = json.dumps(value, sort_keys=True, separators=(",", ":")).encode("utf-8")
    return hashlib.sha256(payload).hexdigest()
