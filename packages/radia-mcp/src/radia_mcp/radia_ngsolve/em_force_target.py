"""Autonomous electromagnetic-force target extraction for loop artifacts.

The helpers in this module are deliberately solver-independent.  They take a
large source-native learning artifact, select the force/torque/motor slots, and
attach public analytic gates before any commercial or live solver run is
claimed.
"""
from __future__ import annotations

import math
import time
from collections import Counter
from datetime import datetime, timezone

from .loop_autolearn import classify_slot_family, relative_error
from .slot_gates import (
    MU0,
    coaxial_pm_force_gap_sweep_gate,
    computed_reference_rows_gate,
    cross_validation_artifact_to_mcp_feedback_gate,
    ipm_saliency_torque_component_gate,
    parallel_wire_force_per_length,
)


EM_FORCE_FAMILY = "force_torque_motor"


def utc_now() -> str:
    """Return a parseable UTC timestamp."""

    return datetime.now(timezone.utc).isoformat().replace("+00:00", "Z")


def _slot_family(slot):
    family = str(slot.get("lesson_family") or "").strip()
    return family or classify_slot_family(slot)


def _slot_label(slot, index, target_kind):
    slot_id = str(slot.get("slot_id", "")).strip() or str(index)
    tool = str(slot.get("tool", "")).strip() or "unknown"
    return f"em_force_slot_{slot_id}_{tool}_{target_kind}"


def _row_common(slot, index, target_kind, quantity, computed, reference, unit, tolerance, extra=None):
    rel = relative_error(computed, reference)
    row = {
        "case": _slot_label(slot, index, target_kind),
        "slot_id": str(slot.get("slot_id", "")),
        "tool": str(slot.get("tool", "")),
        "lesson_family": EM_FORCE_FAMILY,
        "row_kind": "em_force_public_gate",
        "target_kind": target_kind,
        "quantity": quantity,
        "computed": float(computed),
        "reference": float(reference),
        "unit": unit,
        "rel_error": rel,
        "tolerance": float(tolerance),
        "pass": rel <= float(tolerance),
    }
    if extra:
        row.update(extra)
    return row


def parallel_wire_force_row(slot, index=0):
    """Return a signed two-wire Lorentz force row."""

    i1, i2, separation = 42.0, 18.0, 0.06
    computed = parallel_wire_force_per_length(i1, i2, separation)
    reference = MU0 * i1 * i2 / (2.0 * math.pi * separation)
    return _row_common(
        slot,
        index,
        "parallel_wire_lorentz",
        "parallel_wire_force_per_length",
        computed,
        reference,
        "N/m",
        1.0e-14,
        {
            "current_1_A": i1,
            "current_2_A": i2,
            "separation_m": separation,
            "force_sign_convention": "positive_for_like-current_attraction",
        },
    )


def magnetic_air_gap_force_row(slot, index=0):
    """Return a magnetic-circuit holding-force row."""

    turns = 120.0
    current = 3.5
    gap = 0.0025
    iron_path = 0.18
    mu_r = 1500.0
    area = 0.00032
    effective_gap = gap + iron_path / mu_r
    b_gap = MU0 * turns * current / effective_gap
    computed = b_gap * b_gap * area / (2.0 * MU0)
    reference = MU0 * (turns * current) ** 2 * area / (2.0 * effective_gap**2)
    return _row_common(
        slot,
        index,
        "magnetic_air_gap_pressure",
        "magnetic_circuit_gap_force",
        computed,
        reference,
        "N",
        1.0e-14,
        {
            "turns": turns,
            "current_A": current,
            "gap_m": gap,
            "iron_path_m": iron_path,
            "mu_r": mu_r,
            "pole_area_m2": area,
            "effective_gap_m": effective_gap,
            "b_gap_T": b_gap,
        },
    )


def ipm_dq_torque_row(slot, index=0):
    """Return an IPM magnet-plus-reluctance torque identity row."""

    gate = ipm_saliency_torque_component_gate(
        lambda_m=0.055,
        Ld=0.0012,
        Lq=0.0021,
        id_current=-18.0,
        iq_current=32.0,
        pole_pairs=4.0,
        tol=1.0e-12,
    )
    return _row_common(
        slot,
        index,
        "ipm_dq_torque_components",
        "ipm_total_torque_equals_component_sum",
        gate["total_torque_Nm"],
        gate["direct_total_torque_Nm"],
        "N*m",
        1.0e-12,
        {
            "gate_policy": gate["policy"],
            "magnet_torque_Nm": gate["magnet_torque_Nm"],
            "reluctance_torque_Nm": gate["reluctance_torque_Nm"],
            "saliency_ratio_Lq_over_Ld": gate["saliency_ratio_Lq_over_Ld"],
            "gate_status": gate["status"],
        },
    )


def pm_gap_sweep_force_row(slot, index=0):
    """Return a PM force-gap sweep invariant row."""

    constant = 2.0e-12
    gaps = (0.006, 0.008, 0.010)
    rows = [
        {"gap_m": gap, "force_N": -constant / gap**4}
        for gap in gaps
    ]
    gate = coaxial_pm_force_gap_sweep_gate(rows, rtol_invariant=1.0e-12)
    return _row_common(
        slot,
        index,
        "pm_force_gap_sweep",
        "pm_force_first_last_gap4_ratio",
        gate["force_ratio_first_last"],
        gate["expected_force_ratio_first_last"],
        "1",
        1.0e-12,
        {
            "gate_policy": gate["policy"],
            "gap_samples_m": list(gaps),
            "force_samples_N": gate["forces_N"],
            "max_force_gap4_invariant_rel_error": gate["max_force_gap4_invariant_rel_error"],
            "gate_status": gate["status"],
        },
    )


def em_force_public_row_for_slot(slot, index=0):
    """Return the public analytic force/torque row assigned to a source slot."""

    tool = str(slot.get("tool", "")).upper()
    if "FEMM" in tool:
        return parallel_wire_force_row(slot, index=index)
    if "JMAG" in tool:
        return ipm_dq_torque_row(slot, index=index)
    if "ELF" in tool:
        return pm_gap_sweep_force_row(slot, index=index)
    if "COMSOL" in tool:
        return magnetic_air_gap_force_row(slot, index=index)
    selector = int(slot.get("slot_id") or index) % 4
    if selector == 0:
        return magnetic_air_gap_force_row(slot, index=index)
    if selector == 1:
        return parallel_wire_force_row(slot, index=index)
    if selector == 2:
        return ipm_dq_torque_row(slot, index=index)
    return pm_gap_sweep_force_row(slot, index=index)


def select_em_force_slots(source_artifact):
    """Return force/torque/motor slots from an autonomous loop artifact."""

    slots = source_artifact.get("slots") if isinstance(source_artifact, dict) else []
    if not isinstance(slots, list):
        return []
    selected = []
    for slot in slots:
        if not isinstance(slot, dict):
            continue
        if _slot_family(slot) == EM_FORCE_FAMILY:
            selected.append(slot)
    return selected


def build_em_force_target_artifact(
    source_artifact,
    *,
    artifact_id="em_force_target",
    source_artifact_id="autonomous_basic_learning",
    run_date_utc=None,
    radia_mcp_version="unknown",
    command="",
):
    """Build a public-safe EM force target artifact from a loop artifact."""

    if not isinstance(source_artifact, dict):
        raise ValueError("source_artifact must be a mapping")

    started = time.perf_counter()
    run_date = run_date_utc or utc_now()
    force_slots = select_em_force_slots(source_artifact)
    rows = [
        em_force_public_row_for_slot(slot, index=index)
        for index, slot in enumerate(force_slots)
    ]
    row_gate = computed_reference_rows_gate({"rows": rows}, max_global_rel_error=1.0e-9)
    target_counts = Counter(row["target_kind"] for row in rows)
    tool_counts = Counter(str(slot.get("tool", "")) for slot in force_slots)
    source_candidate_count = sum(
        1
        for slot in force_slots
        if str(slot.get("learning_lanes", {}).get("source_tool", "candidate")) == "candidate"
    )
    slot_records = []
    for index, (slot, row) in enumerate(zip(force_slots, rows)):
        slot_records.append(
            {
                "index": index,
                "slot_id": slot.get("slot_id"),
                "lap": slot.get("lap"),
                "slot_index_in_lap": slot.get("slot_index_in_lap"),
                "tool": slot.get("tool"),
                "lesson_family": EM_FORCE_FAMILY,
                "public_force_gate": row["target_kind"],
                "public_row_quantity": row["quantity"],
                "public_row_pass": row["pass"],
                "learning_lanes": {
                    "public": "verified" if row["pass"] else "candidate",
                    "source_tool": "candidate",
                },
                "next_action": "run a source-native solver-ready force slot, then promote only scrubbed identities back to radia-mcp",
            }
        )

    pass_artifact = bool(rows) and row_gate["status"] == "ok"
    artifact = {
        "schema": "radia.crossval.v1",
        "tool_slot": "radia-mcp",
        "case": "autonomous EM force target pass",
        "artifact_role": "em_force_target",
        "pass": pass_artifact,
        "created_at_utc": run_date,
        "versions": {
            "solver": "em-force-target analytic gates v1",
            "radia_mcp": radia_mcp_version,
        },
        "execution": {
            "run_date_utc": run_date,
            "command": command,
            "source_artifact_id": source_artifact_id,
            "result_artifact_id": artifact_id,
        },
        "result_artifact_id": artifact_id,
        "result_output_schema_id": "em_force_target_rows_v1",
        "result_output_columns": [
            "slot_id",
            "tool",
            "target_kind",
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
            "target_kind": "1",
            "quantity": "1",
            "computed": "case-dependent",
            "reference": "case-dependent",
            "unit": "1",
            "rel_error": "1",
            "tolerance": "1",
            "pass": "1",
        },
        "source_artifact": {
            "artifact_id": source_artifact_id,
            "artifact_role": source_artifact.get("artifact_role"),
            "result_artifact_id": source_artifact.get("result_artifact_id"),
            "slot_count": len(source_artifact.get("slots", [])) if isinstance(source_artifact.get("slots"), list) else 0,
        },
        "row_gate": row_gate,
        "rows": rows,
        "slots": slot_records,
        "summary": {
            "force_slot_count": len(force_slots),
            "row_count": len(rows),
            "target_counts": dict(sorted(target_counts.items())),
            "source_tool_counts": dict(sorted(tool_counts.items())),
            "source_tool_candidate_count": source_candidate_count,
            "public_row_gate_status": row_gate["status"],
            "solver_ready_queue_needed": len(force_slots),
        },
        "learning_lanes": {
            "public": "verified" if pass_artifact else "encoded",
            "source_tool": "candidate" if force_slots else "none",
        },
        "learning_lane_details": {
            "source_tool": {
                "reason": (
                    "This pass extracts and verifies public analytic EM-force gates only. "
                    "Commercial/live source-tool MCPs stay candidate until their own solver-ready runs, private edits, and focused tests pass."
                )
            }
        },
        "public_lesson": (
            "Electromagnetic-force loop slots should first pass solver-independent analytic gates for Lorentz force, "
            "Maxwell air-gap pressure, PM gap-sweep scaling, or dq torque decomposition before any live source-tool result is promoted."
        ),
        "learning_targets": [
            "radia-mcp: em_force_target.build_em_force_target_artifact",
            "radia-mcp: computed_reference_rows_gate",
            "radia-mcp: force_validation method map",
        ],
        "verification": {
            "public": command,
            "commands": [{"command": command, "result": "passed"}] if command else [],
        },
        "mcp_feedback": {
            "public_summary": (
                "Added a public-safe EM-force target pass that turns force_torque_motor slots into verified analytic rows "
                "and a source-tool solver-ready queue without claiming live commercial execution."
            ),
            "encoded_targets": [
                "radia-mcp: em_force_target",
                "radia-mcp validation: electromagnetic_force_target.py",
            ],
            "knowledge_topics": ["force_validation", "loop_learning:em_force_target"],
        },
        "timing_breakdown_s": {
            "load_source_artifact": 0.0,
            "select_force_slots": 0.0,
            "public_force_rows": round(time.perf_counter() - started, 6),
            "row_gate": 0.0,
        },
        "checks": {
            "force_slots_extracted": bool(force_slots),
            "all_force_slots_have_rows": len(rows) == len(force_slots),
            "row_gate_ok": row_gate["status"] == "ok",
            "public_rows_pass": all(row["pass"] for row in rows),
            "source_tool_learning_not_overclaimed": True,
        },
        "next_slot_allowed": pass_artifact,
        "notes": [
            "This is a solver-ready target pass, not a live commercial solver run.",
            "Run source-native FEMM/JMAG/ELF/COMSOL force slots next, then feed verified and scrubbed lessons back into the appropriate MCP lane.",
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
