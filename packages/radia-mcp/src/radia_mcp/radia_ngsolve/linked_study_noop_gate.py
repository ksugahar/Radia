"""Solver-neutral evidence gate for silent no-op linked-study runs."""

from __future__ import annotations

import math
from collections.abc import Mapping
from typing import Any


def linked_study_silent_noop_gate(
    summary: Mapping[str, Any],
    *,
    maximum_noop_seconds: float = 1.0,
) -> dict[str, object]:
    """Confirm that a linked-study execution was a measured silent no-op.

    This gate validates a negative result. A native API returning without an
    exception is not solver evidence when every linked study still has no
    result identity, result file, or result table.
    """

    rows = summary.get("studies")
    if not isinstance(rows, list) or len(rows) < 2 or not all(isinstance(row, Mapping) for row in rows):
        raise ValueError("studies must contain at least two study mappings")
    limit = float(maximum_noop_seconds)
    if not math.isfinite(limit) or limit <= 0.0:
        raise ValueError("maximum_noop_seconds must be finite and positive")

    indices = [row.get("index") for row in rows]
    durations: list[float] = []
    try:
        durations = [float(row["run_seconds"]) for row in rows]
    except (KeyError, TypeError, ValueError) as exc:
        raise ValueError("each study must record finite run_seconds") from exc

    digest_before = str(summary.get("source_digest_before") or "")
    digest_after = str(summary.get("source_digest_after") or "")
    checks = {
        "native_hidden_execution_recorded": summary.get("execution_mode") == "native_hidden_linked_study_run",
        "source_digest_preserved": len(digest_before) == 64 and digest_before == digest_after,
        "temporary_work_copy_used": summary.get("work_copy_used") is True,
        "owned_process_released": summary.get("owned_process_released") is True,
        "no_owned_process_left": summary.get("owned_process_count_after") == 0,
        "study_order_contiguous": indices == list(range(len(rows))),
        "study_roles_recorded": all(bool(str(row.get("role") or "").strip()) for row in rows),
        "all_started_without_results": all(row.get("has_result_before") is False for row in rows),
        "all_returned_without_results": all(row.get("has_result_after") is False for row in rows),
        "all_runs_returned_quickly": all(math.isfinite(value) and 0.0 <= value <= limit for value in durations),
        "no_result_files_created": all(int(row.get("result_file_count", -1)) == 0 for row in rows),
        "no_result_tables_created": all(int(row.get("result_table_count", -1)) == 0 for row in rows),
        "solver_success_not_claimed": summary.get("solver_success") is False,
    }
    verified_noop = all(checks.values())
    return {
        "policy": "linked_study_silent_noop_gate_v1",
        "status": "ok" if verified_noop else "needs_attention",
        "classification": "verified_silent_noop" if verified_noop else "incomplete_or_contradictory_evidence",
        "solver_result_accepted": False,
        "checks": checks,
        "issues": [name for name, ok in checks.items() if not ok],
        "metrics": {
            "study_count": len(rows),
            "maximum_run_seconds": max(durations),
            "total_run_seconds": sum(durations),
        },
        "lesson": (
            "Treat API return as execution evidence only. Promote a linked solve to numerical evidence "
            "after every required study has a result identity and the expected result files or tables."
        ),
    }
