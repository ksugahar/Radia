"""Solver-neutral six-degree-of-freedom motion-table validation."""
from __future__ import annotations

import math


def motion_table_coordinate_gate(
    translation_times_s,
    translation_vectors,
    rotation_times_s,
    rotation_vectors,
    *,
    coordinate_frame_id: str,
    translation_unit: str = "mm",
    rotation_unit: str = "deg",
    motion_semantics: str = "cumulative_displacement",
):
    """Validate independent translation and rotation time tables.

    Translation and rotation may legitimately end at different times.  The
    gate therefore validates each axis independently instead of silently
    resampling both onto one synthetic time vector.
    """

    tt = [float(value) for value in translation_times_s]
    tr = [[float(value) for value in row] for row in translation_vectors]
    rt = [float(value) for value in rotation_times_s]
    rr = [[float(value) for value in row] for row in rotation_vectors]
    frame = str(coordinate_frame_id).strip()
    finite = all(math.isfinite(value) for value in tt + rt)
    finite = finite and all(math.isfinite(value) for row in tr + rr for value in row)
    checks = {
        "coordinate_frame_recorded": bool(frame),
        "translation_unit_supported": translation_unit in {"m", "mm"},
        "rotation_unit_supported": rotation_unit in {"rad", "deg"},
        "cumulative_displacement_semantics": motion_semantics == "cumulative_displacement",
        "translation_has_at_least_two_rows": len(tt) >= 2,
        "rotation_has_at_least_two_rows": len(rt) >= 2,
        "translation_row_count_matches": len(tt) == len(tr),
        "rotation_row_count_matches": len(rt) == len(rr),
        "translation_vectors_are_3d": all(len(row) == 3 for row in tr),
        "rotation_vectors_are_3d": all(len(row) == 3 for row in rr),
        "all_finite": finite,
        "translation_time_starts_at_zero": bool(tt) and tt[0] == 0.0,
        "rotation_time_starts_at_zero": bool(rt) and rt[0] == 0.0,
        "translation_time_strictly_increases": all(b > a for a, b in zip(tt, tt[1:])),
        "rotation_time_strictly_increases": all(b > a for a, b in zip(rt, rt[1:])),
    }
    return {
        "policy": "motion_table_coordinate_gate_v1",
        "status": "ok" if all(checks.values()) else "needs_attention",
        "coordinate_frame_id": frame,
        "translation_unit": translation_unit,
        "rotation_unit": rotation_unit,
        "motion_semantics": motion_semantics,
        "translation_sample_count": len(tt),
        "rotation_sample_count": len(rt),
        "translation_end_time_s": tt[-1] if tt else None,
        "rotation_end_time_s": rt[-1] if rt else None,
        "independent_time_axes": tt != rt,
        "combined_motion_end_time_s": max(tt[-1], rt[-1]) if tt and rt else None,
        "checks": checks,
        "lesson": (
            "Keep translation and rotation time axes independent, record the coordinate frame, "
            "units, and cumulative-versus-incremental semantics, and never invent synchronized "
            "rows before interpolation or solver handoff."
        ),
    }
