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
            and fidelity.get("stimulus") == "source-faithful"
            and fidelity.get("boundary_operator") == "source_homogeneous_dirichlet"
            and fidelity.get("geometry") == "source_profile_polygonized"
            and fidelity.get("material") == "source_region_materials"
            and fidelity.get("required_features_mapped") is True
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
    retirement_ready = bool(records) and passed_count == faithful_count == len(records)
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
        "source_solver_not_launched": evidence.get("source_solver_launched") is False,
        "retirement_claim_matches_fidelity": evidence.get("retirement_ready")
        is retirement_ready,
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
