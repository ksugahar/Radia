"""Producer-side validation for nonlinear magnet-model handoffs."""
from __future__ import annotations

import math
import re


_SHA256 = re.compile(r"^(?:sha256:)?[0-9a-fA-F]{64}$")


def magnet_model_handoff_gate(
    residual_phases,
    *,
    nonlinear_tolerance: float,
    source_result_artifact_id: str,
    source_result_digest: str,
    magnet_control_artifact_id: str,
    magnet_control_digest: str,
    magnet_geometry_artifact_id: str,
    magnet_geometry_digest: str,
    numbering_policy: str,
    element_id_offset: int,
    node_id_offset: int,
    material_mapping_count: int,
    geometry_transform: str,
):
    """Check convergence and two-file model identity before downstream use."""

    phases = [[float(value) for value in phase] for phase in residual_phases]
    tolerance = float(nonlinear_tolerance)
    finite = all(math.isfinite(value) and value >= 0.0 for phase in phases for value in phase)
    phase_decrease = all(
        phase and all(right < left for left, right in zip(phase, phase[1:]))
        for phase in phases
    )
    final_residual = phases[-1][-1] if phases and phases[-1] else math.inf
    offsets_are_int = isinstance(element_id_offset, int) and isinstance(node_id_offset, int)
    numbering_ok = numbering_policy in {"preserve", "offset"} and offsets_are_int
    if numbering_policy == "preserve":
        numbering_ok = numbering_ok and element_id_offset == 0 and node_id_offset == 0
    checks = {
        "at_least_one_residual_phase": bool(phases),
        "every_phase_has_samples": bool(phases) and all(bool(phase) for phase in phases),
        "residuals_finite_nonnegative": finite,
        "residual_decreases_within_each_phase": phase_decrease,
        "nonlinear_tolerance_positive": math.isfinite(tolerance) and tolerance > 0.0,
        "terminal_residual_meets_tolerance": final_residual <= tolerance,
        "source_result_artifact_id_recorded": bool(str(source_result_artifact_id).strip()),
        "source_result_digest_is_sha256": bool(_SHA256.fullmatch(str(source_result_digest).strip())),
        "magnet_control_artifact_id_recorded": bool(str(magnet_control_artifact_id).strip()),
        "magnet_control_digest_is_sha256": bool(_SHA256.fullmatch(str(magnet_control_digest).strip())),
        "magnet_geometry_artifact_id_recorded": bool(str(magnet_geometry_artifact_id).strip()),
        "magnet_geometry_digest_is_sha256": bool(_SHA256.fullmatch(str(magnet_geometry_digest).strip())),
        "numbering_policy_and_offsets_consistent": numbering_ok,
        "material_mapping_recorded": isinstance(material_mapping_count, int) and material_mapping_count > 0,
        "geometry_transform_recorded": geometry_transform in {"identity", "translated", "rotated", "translated_and_rotated"},
    }
    return {
        "policy": "magnet_model_handoff_gate_v1",
        "status": "ok" if all(checks.values()) else "needs_attention",
        "handoff_ready": all(checks.values()),
        "phase_count": len(phases),
        "phase_final_residuals": [phase[-1] if phase else None for phase in phases],
        "terminal_residual": final_residual,
        "nonlinear_tolerance": tolerance,
        "source_result_artifact_id": source_result_artifact_id,
        "magnet_control_artifact_id": magnet_control_artifact_id,
        "magnet_geometry_artifact_id": magnet_geometry_artifact_id,
        "numbering_policy": numbering_policy,
        "element_id_offset": element_id_offset,
        "node_id_offset": node_id_offset,
        "material_mapping_count": material_mapping_count,
        "geometry_transform": geometry_transform,
        "checks": checks,
        "lesson": (
            "A downstream permanent-magnet model must bind both control/material and geometry artifacts "
            "to the converged source result. Preserve or explicitly offset element/node identities and "
            "record material remapping plus any geometry transform."
        ),
    }
