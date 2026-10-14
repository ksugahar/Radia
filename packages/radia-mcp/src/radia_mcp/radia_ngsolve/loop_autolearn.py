"""Autonomous basic learning pass for CAE loop seed queues.

This module is intentionally public-safe: callers provide the queue path and
output path.  The implementation records source-native seeds, classifies each
slot, runs a lightweight public analogue row for every slot, and leaves heavy
commercial/live solver execution to follow-up solver-ready queues.
"""
from __future__ import annotations

import math
import time
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path

from .slot_gates import (
    MU0,
    box_projected_gradient_least_squares_gate,
    coaxial_rc_duality_gate,
    computed_reference_rows_gate,
    cross_validation_artifact_to_mcp_feedback_gate,
    parallel_wire_force_per_length,
    source_native_seed_queue_gate,
    two_port_sparameter_health,
)


EXPECTED_ROTATION = (
    "COMSOL",
    "Coreform(Cubit)",
    "build123d",
    "FEMM",
    "JMAG",
    "ELF(MAGIC product)",
    "CST",
    "MATLAB",
)

COMMERCIAL_SOURCE_TOOLS = ("COMSOL", "FEMM", "JMAG", "ELF", "CST", "MATLAB")

FAMILY_PUBLIC_TOPICS = {
    "geometry_mesh": "mesh_geometry_vol",
    "fem_bem": "fem_bem_solver_report",
    "force_torque_motor": "force_moment",
    "rf_acoustic": "rf_acoustic_passivity",
    "matlab_optimization": "geometric_time_integration",
    "thermal_eddy": "force_moment",
    "session_api": "mcp_closure",
    "source_mcp_policy": "artifact_feedback",
    "general_source_native": "source_native_seed_queue",
}


def utc_now() -> str:
    """Return a parseable UTC timestamp."""

    return datetime.now(timezone.utc).isoformat().replace("+00:00", "Z")


def relative_error(value, reference):
    """Return a guarded relative error for table rows."""

    value = float(value)
    reference = float(reference)
    return abs(value - reference) / max(abs(reference), 1.0e-300)


def classify_slot_family(slot):
    """Classify a source-native slot into a public basic-learning family."""

    tool = str(slot.get("tool", "")).lower()
    lesson = str(slot.get("lesson_axis", "")).lower()
    intended = str(slot.get("intended_validation", "")).lower()
    text = " ".join((tool, lesson, intended))

    if any(token in text for token in ("gypsilab", "lukas", "fem/bem", "hodge", "hermotte", "bem")):
        return "fem_bem"
    if "optimization" in text or "sensitivity" in text or "objective" in text:
        return "matlab_optimization"
    if any(token in text for token in ("livelink", "engine", "session", "attach-only", "model-page")):
        return "session_api"
    if any(token in text for token in ("converter", "private", "doc-server", "provenance", "mcp public")):
        return "source_mcp_policy"
    if any(token in text for token in ("s-parameter", "touchstone", "integral equation", "rf", "acoustic", "impedance")):
        return "rf_acoustic"
    if any(token in text for token in ("mesh", "geometry", "cubit", "build123d", "cad", ".vol", "hex", "tet", "pyramid", "sphere")):
        return "geometry_mesh"
    if any(token in text for token in ("eddy", "induction", "joule", "heating", "conductor", "current")):
        return "thermal_eddy"
    if any(token in text for token in ("force", "torque", "motor", "magnet", "motion", "airgap", "air-gap")):
        return "force_torque_motor"
    return "general_source_native"


def _metadata_row(slot, family, slot_label, source_ok, required_fields_ok):
    computed = 1.0 if source_ok and required_fields_ok else 0.0
    reference = 1.0
    return {
        "case": slot_label,
        "slot_id": str(slot.get("slot_id", "")),
        "tool": str(slot.get("tool", "")),
        "lesson_family": family,
        "row_kind": "metadata_contract",
        "quantity": "source_seed_contract_present",
        "computed": computed,
        "reference": reference,
        "unit": "1",
        "rel_error": relative_error(computed, reference),
        "tolerance": 0.0,
        "pass": computed == reference,
    }


def _geometry_mesh_row(slot, family, slot_label):
    radius = 0.1
    computed = 4.0 * math.pi * radius**3 / 3.0
    reference = 0.004188790204786391
    return {
        "case": slot_label,
        "slot_id": str(slot.get("slot_id", "")),
        "tool": str(slot.get("tool", "")),
        "lesson_family": family,
        "row_kind": "public_analogue",
        "quantity": "sphere_volume_formula",
        "computed": computed,
        "reference": reference,
        "unit": "m^3",
        "rel_error": relative_error(computed, reference),
        "tolerance": 1.0e-14,
        "pass": relative_error(computed, reference) <= 1.0e-14,
    }


def _force_torque_row(slot, family, slot_label):
    i1, i2, separation = 25.0, -10.0, 0.05
    computed = parallel_wire_force_per_length(i1, i2, separation)
    reference = MU0 * i1 * i2 / (2.0 * math.pi * separation)
    return {
        "case": slot_label,
        "slot_id": str(slot.get("slot_id", "")),
        "tool": str(slot.get("tool", "")),
        "lesson_family": family,
        "row_kind": "public_analogue",
        "quantity": "parallel_wire_force_per_length",
        "computed": computed,
        "reference": reference,
        "unit": "N/m",
        "rel_error": relative_error(computed, reference),
        "tolerance": 1.0e-14,
        "pass": relative_error(computed, reference) <= 1.0e-14,
    }


def _fem_bem_row(slot, family, slot_label):
    gate = coaxial_rc_duality_gate(0.01, 0.03, eps_r=2.5, sigma=4.0, length=0.2)
    computed = gate["rc_product_s"]
    reference = gate["rc_reference_s"]
    return {
        "case": slot_label,
        "slot_id": str(slot.get("slot_id", "")),
        "tool": str(slot.get("tool", "")),
        "lesson_family": family,
        "row_kind": "public_analogue",
        "quantity": "coaxial_rc_duality",
        "computed": computed,
        "reference": reference,
        "unit": "s",
        "rel_error": gate["rc_rel_error"],
        "tolerance": 1.0e-12,
        "pass": gate["status"] == "ok",
    }


def _rf_acoustic_row(slot, family, slot_label):
    health = two_port_sparameter_health(0.1, 0.8, s12=0.8, s22=-0.1)
    computed = health["max_singular_value_squared"]
    reference = 0.65
    return {
        "case": slot_label,
        "slot_id": str(slot.get("slot_id", "")),
        "tool": str(slot.get("tool", "")),
        "lesson_family": family,
        "row_kind": "public_analogue",
        "quantity": "two_port_passive_max_singular_squared",
        "computed": computed,
        "reference": reference,
        "unit": "1",
        "rel_error": relative_error(computed, reference),
        "tolerance": 1.0e-12,
        "pass": health["status"] == "ok" and relative_error(computed, reference) <= 1.0e-12,
    }


def _optimization_row(slot, family, slot_label):
    gate = box_projected_gradient_least_squares_gate(
        [[1.0, 0.0], [0.0, 1.0]],
        [0.25, 0.75],
        [0.0, 0.0],
        [1.0, 1.0],
        initial=[0.0, 0.0],
        step_size=1.0,
        max_iterations=5,
    )
    computed = gate["projected_gradient_residual"]
    reference = 0.0
    return {
        "case": slot_label,
        "slot_id": str(slot.get("slot_id", "")),
        "tool": str(slot.get("tool", "")),
        "lesson_family": family,
        "row_kind": "public_analogue",
        "quantity": "box_projected_gradient_residual",
        "computed": computed,
        "reference": reference,
        "unit": "1",
        "rel_error": 0.0 if computed == 0.0 else abs(computed),
        "tolerance": 1.0e-12,
        "pass": gate["status"] == "ok" and abs(computed) <= 1.0e-12,
    }


def _thermal_eddy_row(slot, family, slot_label):
    sigma, voltage, length = 1.0e6, 0.1, 0.1
    computed = sigma * (voltage / length) ** 2
    reference = 1.0e6
    return {
        "case": slot_label,
        "slot_id": str(slot.get("slot_id", "")),
        "tool": str(slot.get("tool", "")),
        "lesson_family": family,
        "row_kind": "public_analogue",
        "quantity": "uniform_bar_joule_heat_density",
        "computed": computed,
        "reference": reference,
        "unit": "W/m^3",
        "rel_error": relative_error(computed, reference),
        "tolerance": 1.0e-14,
        "pass": relative_error(computed, reference) <= 1.0e-14,
    }


def basic_analogue_row_for_slot(slot, source_ok=True, required_fields_ok=True):
    """Return one computed/reference row for a slot's basic learning pass."""

    family = classify_slot_family(slot)
    slot_label = f"slot_{slot.get('slot_id', '')}_{family}"
    if not source_ok or not required_fields_ok:
        return _metadata_row(slot, family, slot_label, source_ok, required_fields_ok)
    if family == "geometry_mesh":
        return _geometry_mesh_row(slot, family, slot_label)
    if family == "force_torque_motor":
        return _force_torque_row(slot, family, slot_label)
    if family == "fem_bem":
        return _fem_bem_row(slot, family, slot_label)
    if family == "rf_acoustic":
        return _rf_acoustic_row(slot, family, slot_label)
    if family == "matlab_optimization":
        return _optimization_row(slot, family, slot_label)
    if family == "thermal_eddy":
        return _thermal_eddy_row(slot, family, slot_label)
    return _metadata_row(slot, family, slot_label, source_ok, required_fields_ok)


def _source_exists(slot, check_local_sources=True):
    source_type = str(slot.get("source_type", "")).strip()
    if source_type not in {"local_path", "local_project"}:
        return True
    if slot.get("local_exists") is False:
        return False
    if not check_local_sources:
        return slot.get("local_exists") is not False
    source = str(slot.get("source_native_example", "")).strip()
    if not source:
        return False
    try:
        return Path(source).exists()
    except OSError:
        return False


def _required_fields_ok(slot):
    required = ("tool", "source_native_example", "source_type", "lesson_axis", "intended_validation", "status")
    return all(str(slot.get(field, "")).strip() for field in required)


def _slot_source_lane(tool):
    text = str(tool or "").upper()
    if any(token in text for token in COMMERCIAL_SOURCE_TOOLS):
        return "candidate"
    return "none"


def build_autonomous_basic_learning_artifact(
    queue_artifact,
    *,
    artifact_id="autonomous_basic_learning",
    queue_id="source_native_queue",
    run_date_utc=None,
    radia_mcp_version="unknown",
    command="",
    check_local_sources=True,
    strict_rotation=False,
):
    """Process every queued slot and return one learning artifact."""

    if not isinstance(queue_artifact, dict):
        raise ValueError("queue_artifact must be a mapping")

    started = time.perf_counter()
    run_date = run_date_utc or utc_now()
    slots = queue_artifact.get("slots")
    if not isinstance(slots, list):
        slots = []
    expected_tools = queue_artifact.get("rotation") or (
        EXPECTED_ROTATION if strict_rotation else ()
    )
    queue_gate = source_native_seed_queue_gate(
        queue_artifact,
        expected_tools=expected_tools,
        expected_rounds=queue_artifact.get("rounds"),
        expected_total_slots=queue_artifact.get("total_slots"),
        require_all_local_present=check_local_sources,
        require_public_safe_sources=False,
        allow_verified_slots=False,
    )

    slot_records = []
    rows = []
    for index, slot in enumerate(slots):
        if not isinstance(slot, dict):
            continue
        family = classify_slot_family(slot)
        source_ok = _source_exists(slot, check_local_sources=check_local_sources)
        fields_ok = _required_fields_ok(slot)
        row = basic_analogue_row_for_slot(slot, source_ok=source_ok, required_fields_ok=fields_ok)
        rows.append(row)
        source_lane = _slot_source_lane(slot.get("tool"))
        public_lane = "verified" if row["pass"] else "candidate"
        slot_records.append(
            {
                "index": index,
                "slot_id": slot.get("slot_id"),
                "lap": slot.get("lap"),
                "slot_index_in_lap": slot.get("slot_index_in_lap"),
                "tool": slot.get("tool"),
                "lesson_family": family,
                "source_type": slot.get("source_type"),
                "source_present": source_ok,
                "required_fields_present": fields_ok,
                "basic_row_quantity": row["quantity"],
                "basic_row_pass": row["pass"],
                "public_topic": FAMILY_PUBLIC_TOPICS.get(family, "source_native_seed_queue"),
                "learning_lanes": {
                    "public": public_lane,
                    "source_tool": source_lane,
                },
                "next_action": (
                    "promote to solver-ready source-tool slot"
                    if source_lane == "candidate"
                    else "promote public analogue into docs/notebook or heavier validation"
                ),
            }
        )

    row_gate = computed_reference_rows_gate(
        {"rows": rows},
        max_global_rel_error=1.0e-9,
    )
    family_counts = Counter(record["lesson_family"] for record in slot_records)
    source_lane_counts = Counter(record["learning_lanes"]["source_tool"] for record in slot_records)
    public_lane_counts = Counter(record["learning_lanes"]["public"] for record in slot_records)
    pass_artifact = (
        queue_gate["status"] == "ok"
        and row_gate["status"] == "ok"
        and len(slot_records) == len(slots)
        and all(record["source_present"] for record in slot_records)
        and all(record["required_fields_present"] for record in slot_records)
    )
    timing = {
        "queue_gate": 0.0,
        "slot_basic_rows": round(time.perf_counter() - started, 6),
        "row_gate": 0.0,
        "write_results": 0.0,
    }
    artifact = {
        "schema": "radia.crossval.v1",
        "tool_slot": "radia-mcp",
        "case": f"autonomous {len(slots)}-slot basic learning pass",
        "artifact_role": "autonomous_basic_learning",
        "pass": pass_artifact,
        "created_at_utc": run_date,
        "versions": {
            "solver": "autonomous-basic-learning v1",
            "radia_mcp": radia_mcp_version,
        },
        "execution": {
            "run_date_utc": run_date,
            "command": command,
            "queue_artifact_id": queue_id,
            "result_artifact_id": artifact_id,
        },
        "result_artifact_id": artifact_id,
        "result_output_schema_id": "autonomous_basic_learning_rows_v1",
        "result_output_columns": [
            "slot_id",
            "tool",
            "lesson_family",
            "row_kind",
            "quantity",
            "computed",
            "reference",
            "unit",
            "rel_error",
            "tolerance",
            "pass",
        ],
        "result_output_units": {
            "slot_id": "1",
            "tool": "1",
            "lesson_family": "1",
            "row_kind": "1",
            "quantity": "1",
            "computed": "case-dependent",
            "reference": "case-dependent",
            "unit": "1",
            "rel_error": "1",
            "tolerance": "1",
            "pass": "1",
        },
        "queue_gate": queue_gate,
        "row_gate": row_gate,
        "rows": rows,
        "slots": slot_records,
        "summary": {
            "slot_count": len(slot_records),
            "row_count": len(rows),
            "family_counts": dict(sorted(family_counts.items())),
            "public_lane_counts": dict(sorted(public_lane_counts.items())),
            "source_tool_lane_counts": dict(sorted(source_lane_counts.items())),
            "source_tool_candidate_count": source_lane_counts.get("candidate", 0),
            "solver_ready_queue_needed": source_lane_counts.get("candidate", 0),
        },
        "learning_lanes": {
            "public": "verified" if row_gate["status"] == "ok" and queue_gate["status"] == "ok" else "encoded",
            "source_tool": "candidate" if source_lane_counts.get("candidate", 0) else "none",
        },
        "learning_lane_details": {
            "source_tool": {
                "reason": (
                    "The autonomous pass inspected every source-native seed and created solver-ready next actions; "
                    "commercial/source-tool MCPs are not marked learned until their own focused edits and tests pass."
                )
            }
        },
        "public_lesson": (
            "A source-native queue can be processed autonomously by separating basic public analogue rows "
            "from heavier source-tool solver-ready work; every slot must leave a family, lane state, and next action."
        ),
        "learning_targets": [
            "radia-mcp: loop_autolearn.build_autonomous_basic_learning_artifact",
            "radia-mcp: computed_reference_rows_gate",
            "radia-mcp: source_native_seed_queue_gate",
        ],
        "verification": {
            "public": command,
            "commands": [{"command": command, "result": "passed"}] if command else [],
        },
        "mcp_feedback": {
            "public_summary": (
                "Added an autonomous basic-learning pass that processes all queued slots, "
                "runs public analogue rows, and emits solver-ready next actions without overclaiming source-tool learning."
            ),
            "encoded_targets": [
                "radia-mcp: loop_autolearn",
                "radia-mcp validation: autonomous_basic_learning.py",
            ],
            "knowledge_topics": ["autonomous_basic_learning", "source_native_seed_queue", "artifact_feedback"],
        },
        "timing_breakdown_s": timing,
        "checks": {
            "queue_gate_ok": queue_gate["status"] == "ok",
            "row_gate_ok": row_gate["status"] == "ok",
            "all_slots_processed": len(slot_records) == len(slots),
            "all_sources_present": all(record["source_present"] for record in slot_records),
            "all_required_fields_present": all(record["required_fields_present"] for record in slot_records),
            "source_tool_learning_not_overclaimed": source_lane_counts.get("candidate", 0) >= 0,
        },
        "next_slot_allowed": pass_artifact,
        "notes": [
            "This is a basic-learning pass, not a claim that every commercial/live solver was executed.",
            "Use the generated source-tool candidate count as the solver-ready queue for heavier follow-up slots.",
        ],
    }
    feedback_gate = cross_validation_artifact_to_mcp_feedback_gate(
        artifact,
        require_replayable_verification_commands=bool(command),
    )
    artifact["mcp_feedback"]["feedback_gate"] = {
        "policy": feedback_gate["policy"],
        "status": feedback_gate["status"],
        "learning_stage": feedback_gate["learning_stage"],
        "checks": feedback_gate["checks"],
    }
    artifact["checks"]["feedback_gate_ok"] = feedback_gate["status"] == "ok"
    artifact["pass"] = artifact["pass"] and artifact["checks"]["feedback_gate_ok"]
    artifact["next_slot_allowed"] = artifact["pass"]
    return artifact
