"""Solver-independent two-port sweep artifact summary gate."""
from __future__ import annotations

import json
import math


def rf_sweep_artifact_summary_gate(summary_json: str, passivity_tolerance: float = 1e-3,
                                   reciprocity_tolerance: float = 1e-3) -> dict:
    summary = json.loads(summary_json)
    if not isinstance(summary, dict):
        raise ValueError("summary_json must decode to an object")
    singular = float(summary.get("max_singular_value", math.inf))
    reciprocity = float(summary.get("max_reciprocity_abs", math.inf))
    channels = summary.get("matrix_channels")
    checks = {
        "solver_completed": summary.get("solved") is True,
        "touchstone_is_fresh": summary.get("touchstone_fresh") is True,
        "two_port_extension": summary.get("touchstone_suffix") == ".s2p",
        "frequency_rows_present": isinstance(summary.get("frequency_rows"), int) and summary["frequency_rows"] >= 2,
        "full_two_port_matrix": set(channels or []) == {"S11", "S12", "S21", "S22"},
        "passive": math.isfinite(singular) and singular <= 1.0 + passivity_tolerance,
        "reciprocal": math.isfinite(reciprocity) and reciprocity <= reciprocity_tolerance,
        "solver_version_recorded": bool(str(summary.get("solver_version", "")).strip()),
        "run_date_recorded": bool(str(summary.get("run_date_utc", "")).strip()),
    }
    return {"schema": "radia-rf-sweep-artifact-summary/v1",
            "status": "ok" if all(checks.values()) else "needs_attention",
            "checks": checks, "max_singular_value": singular,
            "max_reciprocity_abs": reciprocity, "frequency_rows": summary.get("frequency_rows")}
