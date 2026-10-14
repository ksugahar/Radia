"""Solver-neutral force-versus-position profile validation."""
from __future__ import annotations

import math


def force_position_profile_gate(
    positions,
    forces,
    *,
    node_counts=None,
    element_counts=None,
    min_sample_count: int = 5,
    max_mesh_count_relative_span: float = 0.02,
    require_interior_peak: bool = False,
    require_nonnegative: bool = False,
):
    x = [float(value) for value in positions]
    f = [float(value) for value in forces]
    if len(x) != len(f):
        raise ValueError("positions and forces must have the same length")
    if min_sample_count < 3:
        raise ValueError("min_sample_count must be >= 3")
    if max_mesh_count_relative_span < 0.0:
        raise ValueError("max_mesh_count_relative_span must be >= 0")

    finite = all(math.isfinite(value) for value in x + f)
    increasing = finite and all(right > left for left, right in zip(x, x[1:]))
    force_range = max(f) - min(f) if f and finite else 0.0
    peak_index = f.index(max(f)) if f and finite else None
    interior_peak = peak_index is not None and 0 < peak_index < len(f) - 1

    def mesh_span(values, name):
        if values is None:
            return None, True
        counts = [int(value) for value in values]
        if len(counts) != len(x):
            raise ValueError(f"{name} must have the same length as positions")
        if not counts or min(counts) <= 0:
            return None, False
        span = (max(counts) - min(counts)) / max(counts)
        return span, span <= max_mesh_count_relative_span

    node_span, node_ok = mesh_span(node_counts, "node_counts")
    element_span, element_ok = mesh_span(element_counts, "element_counts")
    checks = {
        "sample_count_sufficient": len(x) >= min_sample_count,
        "all_finite": finite,
        "positions_strictly_increase": increasing,
        "force_profile_nontrivial": force_range > 0.0,
        "node_count_drift_ok": node_ok,
        "element_count_drift_ok": element_ok,
        "interior_peak_present_when_required": not require_interior_peak or interior_peak,
        "nonnegative_when_required": not require_nonnegative or (finite and min(f) >= 0.0),
    }
    return {
        "policy": "force_position_profile_gate_v1",
        "status": "ok" if all(checks.values()) else "needs_attention",
        "sample_count": len(x),
        "position_min": min(x) if x and finite else None,
        "position_max": max(x) if x and finite else None,
        "force_min": min(f) if f and finite else None,
        "force_max": max(f) if f and finite else None,
        "force_range": force_range,
        "peak_index": peak_index,
        "peak_position": x[peak_index] if peak_index is not None else None,
        "interior_peak": interior_peak,
        "node_count_relative_span": node_span,
        "element_count_relative_span": element_span,
        "max_mesh_count_relative_span": max_mesh_count_relative_span,
        "checks": checks,
        "lesson": (
            "Do not impose monotonicity on an actuator force-position profile. "
            "Preserve its extrema and record remeshing drift alongside the force samples."
        ),
    }
