"""Public-safe helpers for learning from converted-model artifacts.

The CAE loop often revisits models that were already converted from another
tool.  These helpers keep that recovery step honest: a converted metadata stub
is useful teaching material, but it is not a solved mesh; a slot artifact is
closed only when its learning lanes say so explicitly.
"""
from __future__ import annotations

from collections import Counter


VALID_LANE_STATUSES = {"none", "candidate", "encoded", "verified"}


def classify_converted_geometry(geometry_kind):
    """Classify a converted model geometry descriptor.

    ``*_stub`` means the converter preserved metadata and emitted a scaffold,
    but source geometry/mesh still has to be rebuilt or re-exported before the
    model can be treated as a runnable solver case.
    """

    kind = str(geometry_kind or "").strip()
    key = kind.lower()
    is_stub = key.endswith("_stub") or key == "stub" or "stub" in key
    is_vol = key.endswith(".vol") or key in {"netgen_vol", "vol_tri_tet", "tri_tet_vol"}
    return {
        "geometry_kind": kind,
        "is_stub": is_stub,
        "is_tri_tet_vol": is_vol,
        "requires_geometry_rebuild": is_stub,
        "solver_ready": bool(kind) and not is_stub,
    }


def learning_lane_closure(learning_lanes, require_source=False):
    """Return a normalized closure summary for public/source learning lanes."""

    lanes = dict(learning_lanes or {})
    public = lanes.get("public")
    source = lanes.get("source_tool")
    invalid = [
        name for name, value in (("public", public), ("source_tool", source))
        if value is not None and value not in VALID_LANE_STATUSES
    ]
    public_closed = public in {"verified", "none"}
    source_closed = source in {"verified", "none"} or (source == "candidate" and not require_source)
    return {
        "public": public,
        "source_tool": source,
        "invalid_statuses": invalid,
        "public_closed": public_closed,
        "source_closed": source_closed,
        "closed": not invalid and public_closed and source_closed,
        "require_source": bool(require_source),
    }


def converted_model_recovery_summary(records, require_source_slots=()):
    """Summarize converted-model recovery records for a loop dashboard.

    Each record may contain ``tool_slot``, ``pass``, ``geometry_kind``, and
    ``learning_lanes``.  Commercial/private provenance should stay outside
    this public-safe summary.
    """

    required = {str(slot).upper() for slot in require_source_slots}
    rows = []
    for record in records:
        slot = str(record.get("tool_slot") or record.get("tool") or "").strip()
        geom = classify_converted_geometry(record.get("geometry_kind"))
        lanes = learning_lane_closure(
            record.get("learning_lanes"),
            require_source=slot.upper() in required,
        )
        row = {
            "tool_slot": slot,
            "pass": record.get("pass") is True,
            "geometry": geom,
            "learning": lanes,
            "artifact": record.get("artifact") or record.get("path"),
        }
        row["ready_for_reuse"] = row["pass"] and lanes["closed"] and not geom["requires_geometry_rebuild"]
        rows.append(row)

    by_tool = Counter(row["tool_slot"] for row in rows)
    stub_tools = [row["tool_slot"] for row in rows if row["geometry"]["requires_geometry_rebuild"]]
    open_lanes = [
        row["tool_slot"] for row in rows
        if row["pass"] and not row["learning"]["closed"]
    ]
    return {
        "policy": "converted_model_recovery_requires_stub_and_lane_checks",
        "n_records": len(rows),
        "by_tool": dict(by_tool),
        "all_passed": all(row["pass"] for row in rows),
        "stub_geometry_tools": stub_tools,
        "open_learning_lane_tools": open_lanes,
        "ready_for_reuse_count": sum(1 for row in rows if row["ready_for_reuse"]),
        "status": "ok" if rows and all(row["pass"] for row in rows) and not open_lanes else "needs_attention",
        "rows": rows,
    }
