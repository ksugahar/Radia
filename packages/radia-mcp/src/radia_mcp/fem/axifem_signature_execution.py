"""Closed-world gate for axisymmetric H1Henrotte signature executions."""

from __future__ import annotations

import hashlib
import json
import math
import re
from collections.abc import Mapping
from datetime import datetime
from typing import Any


SCHEMA = "radia.axifem-signature-execution.v1"
GATE_SCHEMA = "radia.axifem-signature-execution-gate.v1"
FEATURES = {
    "arc_geometry",
    "circuits",
    "external_region",
    "moving_band",
    "nonlinear_bh",
    "periodic_boundary",
    "point_properties",
}
READINESS = {"source_faithful", "validation_surrogate"}
SHA256_RE = re.compile(r"^[0-9a-f]{64}$")
CASE_RE = re.compile(r"^axisymmetric_dc_case_[0-9]{2}$")
PATH_RE = re.compile(r"(?:^[A-Za-z]:[\\/]|\\|(?:^|/)\.\.(?:/|$))")


def _sha(value: Mapping[str, Any]) -> str:
    payload = dict(value)
    payload.pop("evidence_payload_sha256", None)
    return hashlib.sha256(
        json.dumps(
            payload, ensure_ascii=True, sort_keys=True, separators=(",", ":")
        ).encode("utf-8")
    ).hexdigest()


def _timestamp(value: object) -> bool:
    try:
        parsed = datetime.fromisoformat(str(value).replace("Z", "+00:00"))
    except ValueError:
        return False
    return parsed.tzinfo is not None


def _finite(value: object, *, positive: bool = False, nonnegative: bool = False) -> bool:
    try:
        number = float(value)
    except (TypeError, ValueError):
        return False
    if not math.isfinite(number):
        return False
    if positive:
        return number > 0.0
    if nonnegative:
        return number >= 0.0
    return True


def _feature_map(value: object) -> dict[str, bool] | None:
    if not isinstance(value, Mapping) or set(value) != FEATURES:
        return None
    if not all(isinstance(item, bool) for item in value.values()):
        return None
    return {str(key): bool(item) for key, item in value.items()}


def _constraint_residual_contract(value: object, residual_limit: float) -> bool:
    if not isinstance(value, Mapping):
        return False
    constraint_count = value.get("constraint_count")
    component_count = value.get("constraint_component_count")
    point_count = value.get("point_constraint_count")
    reduced_count = value.get("reduced_dof_count")
    return all(
        (
            isinstance(constraint_count, int) and constraint_count >= 0,
            isinstance(component_count, int) and 0 <= component_count <= constraint_count,
            isinstance(point_count, int) and point_count >= 0,
            isinstance(reduced_count, int) and reduced_count > 0,
            _finite(value.get("max_identification_abs_error"), nonnegative=True)
            and float(value["max_identification_abs_error"]) <= 1.0e-10,
            _finite(value.get("max_point_constraint_abs_error"), nonnegative=True)
            and float(value["max_point_constraint_abs_error"]) <= 1.0e-12,
            _finite(value.get("relative_reduced_residual"), nonnegative=True)
            and float(value["relative_reduced_residual"]) <= residual_limit,
        )
    )


def _point_property_contract(
    source_value: object,
    potential_value: object,
    requested: dict[str, bool] | None,
    mapped: dict[str, bool] | None,
) -> bool:
    if (
        not isinstance(source_value, Mapping)
        or not isinstance(potential_value, Mapping)
        or requested is None
        or mapped is None
    ):
        return False
    source_count = source_value.get("source_count")
    source_embedded = source_value.get("embedded_vertex_count")
    annihilated = source_value.get("axis_annihilated_count")
    potential_count = potential_value.get("constraint_count")
    potential_embedded = potential_value.get("embedded_vertex_count")
    requested_point = requested["point_properties"]
    mapped_point = mapped["point_properties"]
    return all(
        (
            source_value.get("method") == "vertex_dirac_ring_current",
            isinstance(source_count, int) and source_count >= 0,
            source_embedded == source_count,
            isinstance(annihilated, int) and 0 <= annihilated <= source_count,
            _finite(source_value.get("max_vertex_distance_m"), nonnegative=True)
            and float(source_value["max_vertex_distance_m"]) <= 1.0e-10,
            _finite(
                source_value.get("max_weak_load_identity_abs_error_a_m"),
                nonnegative=True,
            )
            and float(source_value["max_weak_load_identity_abs_error_a_m"])
            <= 1.0e-12,
            source_value.get("all_point_properties_exact") is mapped_point,
            potential_value.get("method") == "vertex_essential_aphi",
            isinstance(potential_count, int) and potential_count >= 0,
            potential_embedded == potential_count,
            _finite(potential_value.get("max_vertex_distance_m"), nonnegative=True)
            and float(potential_value["max_vertex_distance_m"]) <= 1.0e-10,
            _finite(
                potential_value.get("max_constraint_abs_error_wb_per_m"),
                nonnegative=True,
            )
            and float(potential_value["max_constraint_abs_error_wb_per_m"])
            <= 1.0e-12,
            (not requested_point and source_count + potential_count == 0)
            or (requested_point and source_count + potential_count > 0),
        )
    )


def _bh_curve_contract(
    value: object,
    nonlinear: bool,
    requested: dict[str, bool] | None,
    mapped: dict[str, bool] | None,
) -> bool:
    if requested is None or mapped is None:
        return False
    if not nonlinear:
        return (
            value in (None, [])
            and requested["nonlinear_bh"] is False
            and mapped["nonlinear_bh"] is False
        )
    if (
        requested["nonlinear_bh"] is not True
        or mapped["nonlinear_bh"] is not True
        or not isinstance(value, list)
        or not value
    ):
        return False
    for row in value:
        if not isinstance(row, Mapping):
            return False
        source_count = row.get("source_point_count")
        if not all(
            (
                bool(SHA256_RE.fullmatch(str(row.get("material_identity_sha256", "")))),
                bool(SHA256_RE.fullmatch(str(row.get("source_curve_sha256", "")))),
                isinstance(source_count, int) and source_count >= 2,
                row.get("effective_point_count") == source_count,
                isinstance(row.get("smoothing_passes"), int)
                and row["smoothing_passes"] >= 0,
                row.get("algorithm") == "source-compatible-natural-cubic-hermite",
                _finite(
                    row.get("max_effective_knot_abs_error_a_per_m"),
                    nonnegative=True,
                )
                and float(row["max_effective_knot_abs_error_a_per_m"]) <= 1.0e-9,
            )
        ):
            return False
    return True


def _boundary_operator_contract(
    value: object,
    requested: dict[str, bool] | None,
    mapped: dict[str, bool] | None,
    fidelity: Mapping[str, Any],
    residual_limit: float,
) -> bool:
    if (
        not isinstance(value, Mapping)
        or requested is None
        or mapped is None
        or value.get("mapped") is not True
    ):
        return False
    operator = str(fidelity.get("boundary_operator", ""))
    if not operator.startswith("source_"):
        return False
    tokens = operator.removeprefix("source_").split("+")
    allowed = {
        "natural",
        "homogeneous_dirichlet",
        "natural_neumann",
        "mixed_robin",
        "dual_boundary_average",
        "signed_periodic_trace",
    }
    if not tokens or any(token not in allowed for token in tokens):
        return False

    source_count = value.get("source_boundary_object_count")
    mixed = value.get("mixed_boundaries")
    periodic = value.get("periodic_pairs")
    residual = value.get("constraint_residual")
    if not all(
        (
            isinstance(source_count, int) and source_count >= 0,
            isinstance(mixed, list),
            isinstance(periodic, list),
        )
    ):
        return False
    for row in mixed:
        if not isinstance(row, Mapping) or not all(
            (
                str(row.get("boundary", "")).startswith("source_mixed_"),
                _finite(row.get("c0")),
                _finite(row.get("c1")),
            )
        ):
            return False
    has_neumann = "natural_neumann" in tokens
    has_robin = "mixed_robin" in tokens
    if (has_neumann or has_robin) != bool(mixed):
        return False
    if has_neumann and any(float(row["c0"]) or float(row["c1"]) for row in mixed):
        return False
    if has_robin and not any(float(row["c0"]) or float(row["c1"]) for row in mixed):
        return False

    periodic_constraints = 0
    for row in periodic:
        if not isinstance(row, Mapping):
            return False
        vertex_count = row.get("vertex_pair_count")
        edge_count = row.get("edge_pair_count")
        constraint_count = row.get("constraint_count")
        if not all(
            (
                isinstance(row.get("boundary_property_index"), int),
                row.get("trace_kind") in {"segment", "arc"},
                row.get("phase") in {-1.0, 1.0},
                row.get("source_object_count") == 2,
                isinstance(vertex_count, int) and vertex_count >= 2,
                isinstance(edge_count, int) and edge_count >= 1,
                isinstance(constraint_count, int) and constraint_count > 0,
                constraint_count <= vertex_count + edge_count,
                _finite(row.get("coordinate_tolerance_m"), nonnegative=True)
                and float(row["coordinate_tolerance_m"]) <= 1.0e-6,
            )
        ):
            return False
        periodic_constraints += constraint_count

    has_periodic = "signed_periodic_trace" in tokens
    if requested["periodic_boundary"] is not has_periodic:
        return False
    if mapped["periodic_boundary"] is not has_periodic:
        return False
    if has_periodic != bool(periodic):
        return False

    has_dual = "dual_boundary_average" in tokens
    if has_dual:
        if value.get("dual_boundary") != "source_dual" or not isinstance(
            residual, Mapping
        ):
            return False
        if not all(
            (
                residual.get("method") == "dual_component_residuals",
                _constraint_residual_contract(residual.get("natural"), residual_limit),
                _constraint_residual_contract(residual.get("essential"), residual_limit),
                _finite(
                    residual.get("average_identity_max_abs_error"),
                    nonnegative=True,
                )
                and float(residual["average_identity_max_abs_error"]) <= 1.0e-12,
            )
        ):
            return False
        active_residual = residual["natural"]
    else:
        if value.get("dual_boundary") is not None or not _constraint_residual_contract(
            residual, residual_limit
        ):
            return False
        active_residual = residual
    return active_residual.get("constraint_count") == periodic_constraints


def _external_region_contract(
    value: object,
    requested: dict[str, bool] | None,
    mapped: dict[str, bool] | None,
) -> bool:
    if (
        not isinstance(value, Mapping)
        or requested is None
        or mapped is None
        or value.get("mapped") is not True
    ):
        return False
    region_count = value.get("region_count")
    parameters = value.get("parameters_m")
    samples = value.get("sample_factors")
    if not all(
        (
            isinstance(region_count, int) and region_count >= 0,
            isinstance(parameters, Mapping),
            isinstance(samples, list),
            not requested["external_region"] or mapped["external_region"] is True,
        )
    ):
        return False
    if not requested["external_region"]:
        return region_count == 0 and not parameters and not samples
    return all(
        (
            region_count > 0,
            set(parameters) == {"z0", "outer_radius", "inner_radius"},
            _finite(parameters.get("z0")),
            _finite(parameters.get("outer_radius"), positive=True),
            _finite(parameters.get("inner_radius"), positive=True),
            bool(samples) and all(_finite(item, positive=True) for item in samples),
            value.get("coefficient_identity")
            == "nu_external=nu_material*(r2+(z-z0)2)*Ri/Ro3",
        )
    )


def _contains_path(value: object, key: str = "") -> bool:
    if key.lower().endswith(("path", "paths", "relative_path")):
        return True
    if isinstance(value, Mapping):
        return any(_contains_path(item, str(name)) for name, item in value.items())
    if isinstance(value, list):
        return any(_contains_path(item, key) for item in value)
    return isinstance(value, str) and bool(PATH_RE.search(value))


def validate_axifem_signature_execution(evidence: Mapping[str, Any]) -> dict[str, Any]:
    """Accept numerical executions while keeping replacement readiness stricter."""

    if not isinstance(evidence, Mapping):
        raise ValueError("evidence must be an object")
    records = evidence.get("records")
    records = records if isinstance(records, list) else []
    failures: list[str] = []
    passed_count = 0
    faithful_count = 0
    surrogate_count = 0
    case_ids: set[str] = set()
    for index, row in enumerate(records):
        label = f"record[{index}]"
        if not isinstance(row, Mapping):
            failures.append(f"{label} is not an object")
            continue
        case_id = str(row.get("case_id", ""))
        requested = _feature_map(row.get("requested_features"))
        mapped = _feature_map(row.get("mapped_features"))
        mesh = row.get("mesh") if isinstance(row.get("mesh"), Mapping) else {}
        display = (
            row.get("display_mesh")
            if isinstance(row.get("display_mesh"), Mapping)
            else {}
        )
        result = row.get("result") if isinstance(row.get("result"), Mapping) else {}
        fidelity = (
            row.get("source_model_fidelity")
            if isinstance(row.get("source_model_fidelity"), Mapping)
            else {}
        )
        point_source = row.get("point_source_evidence")
        point_potential = row.get("point_potential_evidence")
        bh_curves = row.get("bh_curve_evidence")
        boundary_evidence = row.get("boundary_operator_evidence")
        external_evidence = row.get("external_region_evidence")
        timing = (
            row.get("timing_seconds")
            if isinstance(row.get("timing_seconds"), Mapping)
            else {}
        )
        nonlinear = result.get("nonlinear_executed") is True
        residual_limit = 1.0e-4 if nonlinear else 1.0e-7
        required_mapped = bool(
            requested is not None
            and mapped is not None
            and all(not enabled or mapped[name] for name, enabled in requested.items())
        )
        faithful = (
            required_mapped
            and fidelity.get("stimulus")
            in {
                "source-faithful",
                "exact-vertex-dirac-ring-current",
                "exact-vertex-essential-aphi",
                "exact-vertex-current-and-potential",
            }
            and _boundary_operator_contract(
                boundary_evidence,
                requested,
                mapped,
                fidelity,
                residual_limit,
            )
            and _external_region_contract(external_evidence, requested, mapped)
            and fidelity.get("geometry") == "source_profile_polygonized"
            and fidelity.get("material") == "source_region_materials"
            and fidelity.get("required_features_mapped") is True
            and row.get("migration_disposition") == "migrate_now"
            and row.get("migration_blockers") == []
        )
        readiness = str(row.get("readiness_class", ""))
        row_ok = all(
            (
                bool(CASE_RE.fullmatch(case_id)),
                case_id not in case_ids,
                bool(SHA256_RE.fullmatch(str(row.get("private_record_sha256", "")))),
                row.get("solver_lane") == "axifem_henrotte",
                row.get("formulation") == "axisymmetric_magnetostatic_Aphi",
                requested is not None,
                mapped is not None,
                _point_property_contract(
                    point_source, point_potential, requested, mapped
                ),
                _bh_curve_contract(bh_curves, nonlinear, requested, mapped),
                _boundary_operator_contract(
                    boundary_evidence,
                    requested,
                    mapped,
                    fidelity,
                    residual_limit,
                ),
                _external_region_contract(external_evidence, requested, mapped),
                mesh.get("format") == "netgen_vol",
                bool(SHA256_RE.fullmatch(str(mesh.get("sha256", "")))),
                isinstance(mesh.get("element_count"), int)
                and mesh["element_count"] > 0,
                isinstance(mesh.get("vertex_count"), int) and mesh["vertex_count"] > 0,
                isinstance(mesh.get("polygon_count"), int) and mesh["polygon_count"] > 0,
                display.get("format") == "gmsh_msh_4_1",
                bool(SHA256_RE.fullmatch(str(display.get("sha256", "")))),
                result.get("execution_status") == "passed",
                result.get("finite_nonzero_solution") is True,
                _finite(result.get("solution_norm"), positive=True),
                _finite(result.get("field_l2_sq_t2_m3"), positive=True),
                _finite(result.get("relative_algebraic_residual"), nonnegative=True)
                and float(result["relative_algebraic_residual"]) <= residual_limit,
                _finite(result.get("residual_limit"), positive=True)
                and float(result["residual_limit"]) == residual_limit,
                readiness in READINESS,
                (readiness == "source_faithful") == faithful,
                set(timing)
                == {
                    "parse_and_semantics",
                    "mesh_and_vol_roundtrip",
                    "solve_and_verify",
                    "total",
                },
                all(_finite(value, nonnegative=True) for value in timing.values()),
                float(timing.get("total", -1.0))
                >= sum(
                    float(timing.get(name, 0.0))
                    for name in (
                        "parse_and_semantics",
                        "mesh_and_vol_roundtrip",
                        "solve_and_verify",
                    )
                )
                - 1.0e-6,
            )
        )
        case_ids.add(case_id)
        if not row_ok:
            failures.append(f"{label} violates the execution contract")
        else:
            passed_count += 1
            faithful_count += int(faithful)
            surrogate_count += int(not faithful)

    case_count = evidence.get("case_count")
    aggregate_ok = all(
        (
            isinstance(case_count, int) and case_count > 0,
            case_count == len(records),
            evidence.get("execution_passed_count") == passed_count,
            evidence.get("source_faithful_solver_ready_count") == faithful_count,
            evidence.get("validation_surrogate_count") == surrogate_count,
        )
    )
    retirement_ready = (
        len(records) == 27 and passed_count == faithful_count == len(records)
    )
    checks = {
        "schema": evidence.get("schema") == SCHEMA,
        "execution_timestamp_recorded": _timestamp(evidence.get("executed_at_utc")),
        "version_identity_complete": isinstance(evidence.get("execution_version"), Mapping)
        and all(
            str(evidence["execution_version"].get(name, "")).strip()
            for name in ("producer", "producer_version", "radia_mcp", "ngsolve", "gmsh")
        ),
        "public_boundary_has_no_paths": not _contains_path(evidence),
        "record_contracts": not failures,
        "aggregate_counts": aggregate_ok,
        "all_27_signature_representatives_present": len(records) == 27,
        "source_solver_not_launched": evidence.get("source_solver_launched") is False,
        "retirement_claim_matches_fidelity": evidence.get("retirement_ready")
        is retirement_ready,
        "classification_complete": evidence.get("classification_complete") is True,
        "evidence_payload_sha256": str(evidence.get("evidence_payload_sha256", ""))
        == _sha(evidence),
    }
    passed = all(checks.values())
    return {
        "schema": GATE_SCHEMA,
        "status": "accepted" if passed else "rejected",
        "pass": passed,
        "checks": checks,
        "case_count": len(records),
        "execution_passed_count": passed_count,
        "source_faithful_solver_ready_count": faithful_count,
        "validation_surrogate_count": surrogate_count,
        "retirement_ready": retirement_ready if passed else False,
        "record_failures": failures[:20],
    }
