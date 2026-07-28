"""Closed-world evidence gate for the six shipping axifem element paths."""

from __future__ import annotations

import hashlib
import json
import math
from collections.abc import Mapping
from typing import Any


SCHEMA = "radia.axifem-element-evidence.v2"
GATE_SCHEMA = "radia.axifem-element-evidence-gate.v1"
FAMILIES = {
    "P1": ("triangle", 1, False),
    "Q1": ("quadrilateral", 1, False),
    "P2": ("triangle", 2, False),
    "Q2": ("quadrilateral", 2, False),
    "P2_curved": ("triangle", 2, True),
    "Q2_curved": ("quadrilateral", 2, True),
}


def _sha(value: Any) -> str:
    return hashlib.sha256(
        json.dumps(value, sort_keys=True, separators=(",", ":")).encode("utf-8")
    ).hexdigest()


def _digest(value: Any) -> bool:
    text = str(value)
    return len(text) == 64 and all(char in "0123456789abcdef" for char in text)


def _git_sha(value: Any) -> bool:
    text = str(value)
    return len(text) == 40 and all(char in "0123456789abcdef" for char in text)


def _finite(value: Any) -> float:
    result = float(value)
    if not math.isfinite(result):
        raise ValueError("numeric evidence must be finite")
    return result


def validate_axifem_element_evidence(evidence: Mapping[str, Any]) -> dict[str, Any]:
    """Validate current, content-addressed evidence for all axifem paths."""

    if not isinstance(evidence, Mapping):
        raise ValueError("evidence must be an object")
    rows = evidence.get("elements")
    if not isinstance(rows, list):
        raise ValueError("elements must be an array")
    by_family = {
        str(row.get("family")): row for row in rows if isinstance(row, Mapping)
    }
    checks: dict[str, bool] = {
        "schema": evidence.get("schema") == SCHEMA,
        "execution_timestamp_recorded": bool(evidence.get("executed_at_utc")),
        "version_identity_complete": all(
            str(evidence.get("execution_version", {}).get(key, "")).strip()
            for key in ("radia", "ngsolve", "python")
        ),
        "git_head_is_full_sha": _git_sha(evidence.get("git_head")),
        "worktree_clean": evidence.get("git_dirty") is False,
        "vol_route": evidence.get("mesh_route") == "Netgen .vol -> ngsolve.Mesh(path)",
        "test_suite_passed": (
            isinstance(evidence.get("test_summary"), Mapping)
            and int(evidence["test_summary"].get("passed", 0)) >= 45
            and int(evidence["test_summary"].get("failed", -1)) == 0
        ),
        "families_exact": set(by_family) == set(FAMILIES) and len(rows) == len(FAMILIES),
    }
    mesh_hashes: dict[str, str] = {}
    for family, (cell_type, order, curved) in FAMILIES.items():
        row = by_family.get(family, {})
        identity = row.get("identity_error_l2_sq", {}) if isinstance(row, Mapping) else {}
        checks[f"{family}_contract"] = (
            isinstance(row, Mapping)
            and row.get("cell_type") == cell_type
            and row.get("order") == order
            and row.get("curved_geometry") is curved
            and row.get("vol_roundtrip") is True
            and _digest(row.get("vol_sha256"))
            and _digest(row.get("mesh_contract_sha256"))
        )
        try:
            identity_values = [_finite(identity.get(name)) for name in ("interpolation", "gradient", "field")]
        except (TypeError, ValueError):
            identity_values = [math.inf]
        checks[f"{family}_limited_de_rham_identity"] = max(identity_values) <= 1.0e-14
        if isinstance(row, Mapping) and _digest(row.get("mesh_contract_sha256")):
            mesh_hashes[family] = str(row["mesh_contract_sha256"])

    p2 = evidence.get("p2_curved_metrics", {})
    q2 = evidence.get("q2_curved_metrics", {})
    try:
        p2_geometry_ratio = _finite(p2.get("curved_volume_error_percent")) / _finite(
            p2.get("straight_volume_error_percent")
        )
        p2_flux_ratio = _finite(p2.get("curved_total_flux_error_percent")) / _finite(
            p2.get("straight_total_flux_error_percent")
        )
        q2_equivalence = _finite(q2.get("maximum_straight_equivalence_relative_error"))
        changes = [_finite(value) for value in q2.get("annular_successive_changes", [])]
    except (TypeError, ValueError, ZeroDivisionError):
        p2_geometry_ratio = p2_flux_ratio = q2_equivalence = math.inf
        changes = []
    checks["p2_curved_geometry_improves_tenfold"] = abs(p2_geometry_ratio) <= 0.1
    checks["p2_curved_flux_improves_twofold"] = abs(p2_flux_ratio) <= 0.5
    checks["q2_curved_matches_straight"] = q2_equivalence <= 1.0e-9
    checks["q2_curved_annular_converges"] = len(changes) >= 3 and changes[-1] < changes[0]

    payload = dict(evidence)
    declared_sha = payload.pop("evidence_payload_sha256", None)
    calculated_sha = _sha(payload)
    checks["evidence_payload_sha256"] = declared_sha == calculated_sha
    passed = all(checks.values())
    return {
        "schema": GATE_SCHEMA,
        "status": "accepted" if passed else "rejected",
        "pass": passed,
        "checks": checks,
        "evidence_payload_sha256": calculated_sha,
        "element_families": sorted(by_family),
        "mesh_contract_sha256_by_family": dict(sorted(mesh_hashes.items())),
        "limited_de_rham_scope": "Henrotte scalar interpolation/gradient/axisymmetric-field identities; not a full H1-HCurl-HDiv-L2 complex",
    }
