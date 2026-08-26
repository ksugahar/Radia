"""Reusable gates for Cubit Nastran mesh-consumer evidence."""

from __future__ import annotations

import math
import re
from collections.abc import Mapping
from typing import Any


_SHA256 = re.compile(r"^[0-9a-fA-F]{64}$")
_EXPORT_VERB = re.compile(r"^\s*export\s+([a-zA-Z0-9_]+)\b")


def _mapping(value: Any, name: str) -> Mapping[str, Any]:
    if not isinstance(value, Mapping):
        raise ValueError(f"{name} must be an object")
    return value


def _positive_int(value: Any, name: str) -> int:
    if isinstance(value, bool):
        raise ValueError(f"{name} must be a positive integer")
    try:
        number = int(value)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"{name} must be a positive integer") from exc
    if number <= 0:
        raise ValueError(f"{name} must be a positive integer")
    return number


def _finite_nonnegative(value: Any, name: str) -> float:
    try:
        number = float(value)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"{name} must be a finite non-negative number") from exc
    if not math.isfinite(number) or number < 0:
        raise ValueError(f"{name} must be a finite non-negative number")
    return number


def evaluate_nastran_consumer_contract(
    summary: Mapping[str, Any],
    *,
    bbox_tolerance: float = 1e-10,
    require_material_assignment: bool = True,
    require_set_semantics: bool = False,
) -> dict[str, Any]:
    """Separate BDF producer validity from downstream consumer behavior.

    The contract deliberately treats Cubit's BDF as a mesh-interchange
    artifact, not a complete analysis deck.  A downstream importer may parse a
    valid file and still lose order, remesh it, or fail to retain the imported
    mesh in its analysis study; those are consumer limitations, not exporter
    failures.
    """

    summary = _mapping(summary, "summary")
    producer = _mapping(summary.get("producer"), "producer")
    consumer = _mapping(summary.get("consumer"), "consumer")
    tolerance = _finite_nonnegative(bbox_tolerance, "bbox_tolerance")

    command = str(producer.get("export_command", "")).strip()
    match = _EXPORT_VERB.match(command)
    verb = match.group(1).lower() if match else ""
    scope = str(producer.get("artifact_scope", "")).strip().lower()
    digest = str(producer.get("bdf_sha256", "")).strip()
    dimension = _positive_int(producer.get("dimension"), "producer.dimension")
    order = _positive_int(producer.get("order"), "producer.order")
    expected_nodes = _positive_int(
        producer.get("expected_nodes"), "producer.expected_nodes"
    )
    linear_nodes = _positive_int(
        producer.get("expected_linear_nodes", expected_nodes),
        "producer.expected_linear_nodes",
    )
    expected_elements = _positive_int(
        producer.get("expected_primary_elements"),
        "producer.expected_primary_elements",
    )
    if dimension not in (2, 3):
        raise ValueError("producer.dimension must be 2 or 3")
    if order not in (1, 2):
        raise ValueError("producer.order must be 1 or 2")
    if linear_nodes > expected_nodes:
        raise ValueError("producer.expected_linear_nodes must not exceed expected_nodes")

    checks = {
        "canonical_export_command": verb == "nastran_bdf",
        "mesh_interchange_scope": scope == "mesh_interchange",
        "independent_parse_ok": producer.get("independent_parse_ok") is True,
        "digest_recorded": bool(_SHA256.fullmatch(digest)),
        "import_succeeded": consumer.get("import_succeeded") is True,
        "not_remeshed": consumer.get("remeshed") is False,
        "study_retains_imported_mesh": consumer.get("study_has_imported_mesh") is True,
        "dimension_preserved": consumer.get("imported_dimension") == dimension,
        "node_count_preserved": consumer.get("imported_nodes") == expected_nodes,
        "primary_element_count_preserved": (
            consumer.get("imported_primary_elements") == expected_elements
        ),
        "material_assignment_available": (
            not require_material_assignment
            or consumer.get("material_assignment_available") is True
        ),
        "set_semantics_verified": (
            not require_set_semantics or consumer.get("set_semantics_verified") is True
        ),
    }

    has_second_order = consumer.get("has_second_order_elements")
    checks["order_preserved"] = (
        order == 1
        or (
            consumer.get("imported_nodes") == expected_nodes
            and has_second_order is not False
        )
    )

    bbox_error = consumer.get("bbox_max_abs_error")
    if bbox_error is None:
        checks["bbox_verified"] = False
    else:
        checks["bbox_verified"] = (
            _finite_nonnegative(bbox_error, "consumer.bbox_max_abs_error")
            <= tolerance
        )

    warnings: list[str] = []
    if verb == "jmag_nastran":
        warnings.append(
            "export jmag_nastran is a deprecated compatibility alias; "
            "use export nastran_bdf for new journals"
        )
    if order == 2 and has_second_order is None and checks["order_preserved"]:
        warnings.append(
            "second-order preservation is inferred from node count because the "
            "consumer did not expose element-order metadata"
        )

    if verb == "nastran":
        status = "wrong_exporter"
        recommendation = "Use export nastran_bdf; Cubit's built-in export nastran is a different contract."
    elif verb not in ("nastran_bdf", "jmag_nastran"):
        status = "wrong_exporter"
        recommendation = "Record the exact export nastran_bdf command used to create the BDF."
    elif verb == "jmag_nastran":
        status = "legacy_alias"
        recommendation = "Migrate the journal to export nastran_bdf before promotion."
    elif scope != "mesh_interchange":
        status = "scope_mismatch"
        recommendation = "Treat this BDF as mesh interchange, not a complete analysis deck."
    elif not checks["independent_parse_ok"] or not checks["digest_recorded"]:
        status = "producer_invalid"
        recommendation = "Re-run an independent BDF parser and record the consumed file digest."
    elif not checks["import_succeeded"]:
        status = "consumer_import_failed"
        recommendation = "Diagnose the downstream importer without blaming a producer that passed independent parsing."
    elif not checks["not_remeshed"]:
        status = "consumer_remeshed"
        recommendation = "Measure imported counts before any consumer-side CreateMesh/remesh operation."
    elif not checks["study_retains_imported_mesh"]:
        status = "consumer_mesh_not_retained"
        recommendation = "Route this dimension/order through a consumer path that retains imported analysis meshes."
    elif not checks["dimension_preserved"]:
        status = "consumer_dimension_mismatch"
        recommendation = "Reject the import because the consumer changed mesh dimension."
    elif not checks["primary_element_count_preserved"]:
        status = "consumer_element_mismatch"
        recommendation = "Compare primary element counts before material or solver setup."
    elif not checks["order_preserved"]:
        status = "consumer_order_downgrade"
        recommendation = "Do not call the P2 handoff validated; use P1 or another consumer route."
    elif not checks["node_count_preserved"]:
        status = "consumer_node_mismatch"
        recommendation = "Reject the import because the consumed node identity changed."
    elif not checks["bbox_verified"]:
        status = "consumer_bbox_unverified"
        recommendation = "Compare the imported bounding box with independent BDF coordinates."
    elif not checks["material_assignment_available"]:
        status = "consumer_material_assignment_unverified"
        recommendation = "Verify that materials can be assigned after mesh import."
    elif not checks["set_semantics_verified"]:
        status = "consumer_set_semantics_unverified"
        recommendation = "Verify requested set semantics separately from BDF parsing."
    else:
        status = "pass"
        recommendation = "Promote this exact dimension/order/consumer contract with its digest."

    return {
        "policy": "cubit_nastran_consumer_contract_v1",
        "status": status,
        "passed": status == "pass",
        "producer_valid": all(
            checks[name]
            for name in (
                "canonical_export_command",
                "mesh_interchange_scope",
                "independent_parse_ok",
                "digest_recorded",
            )
        ),
        "consumer_valid": status == "pass",
        "dimension": dimension,
        "order": order,
        "expected_nodes": expected_nodes,
        "expected_linear_nodes": linear_nodes,
        "expected_primary_elements": expected_elements,
        "checks": checks,
        "warnings": warnings,
        "recommendation": recommendation,
    }

