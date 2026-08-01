"""Evidence gate for converged HDiv-MMM motor torque on structured HEX meshes."""

from __future__ import annotations

import math
from typing import Any, Mapping, Sequence


_SPACE = {
    "family": "HDiv",
    "order": 1,
    "cell_family": "HEX",
    "project_lane": "BDM1",
    "strict_name": "tensor_product_hdiv_order1",
    "simplex_analogue": "BDM1",
}


def _finite_vector(values: Any, name: str) -> list[float]:
    if not isinstance(values, Sequence) or isinstance(values, (str, bytes)):
        raise ValueError(f"{name} must be a sequence")
    parsed = [float(value) for value in values]
    if len(parsed) < 5 or not all(math.isfinite(value) for value in parsed):
        raise ValueError(f"{name} must contain at least five finite values")
    return parsed


def _relative_l2(left: Sequence[float], right: Sequence[float]) -> float:
    numerator = math.sqrt(sum((a - b) ** 2 for a, b in zip(left, right)))
    denominator = max(math.sqrt(sum(value * value for value in left)), 1.0e-30)
    return numerator / denominator


def hdiv_hex_motor_torque_gate(
    summary: Mapping[str, Any],
    *,
    refinement_tolerance: float = 1.0e-2,
    route_tolerance: float = 2.0e-5,
    reference_tolerance: float = 3.0e-2,
    endpoint_tolerance: float = 1.0e-3,
) -> dict[str, Any]:
    """Validate BDM1-lane/HEX torque convergence and an independent reference."""
    if not isinstance(summary, Mapping):
        raise ValueError("summary must be a mapping")
    levels = summary.get("levels")
    if not isinstance(levels, list) or len(levels) < 3:
        raise ValueError("levels must contain at least three refinement levels")

    parsed_levels = []
    for index, level in enumerate(levels):
        if not isinstance(level, Mapping):
            raise ValueError(f"levels[{index}] must be a mapping")
        moment = _finite_vector(level.get("torque_moment_Nm"), f"levels[{index}].torque_moment_Nm")
        virtual = _finite_vector(
            level.get("torque_virtual_work_Nm"),
            f"levels[{index}].torque_virtual_work_Nm",
        )
        if len(moment) != len(virtual):
            raise ValueError("moment and virtual-work torque arrays must have equal length")
        parsed_levels.append(
            {
                "hex_element_count": int(level.get("hex_element_count", 0)),
                "ndof": int(level.get("ndof", 0)),
                "moment": moment,
                "virtual": virtual,
                "space": dict(level.get("discrete_space") or {}),
                "operator_build_count": int(level.get("operator_build_count", 0)),
            }
        )

    final = parsed_levels[-1]
    prior = parsed_levels[-2]
    if len(final["moment"]) != len(prior["moment"]):
        raise ValueError("the final two refinement levels must share an angle grid")
    reference = summary.get("reference")
    if not isinstance(reference, Mapping):
        raise ValueError("reference must be an executed same-identity mapping")
    reference_torque = _finite_vector(reference.get("torque_Nm"), "reference.torque_Nm")
    if len(reference_torque) != len(final["moment"]):
        raise ValueError("reference and final torque arrays must have equal length")

    peak = max(max(abs(value) for value in final["moment"]), 1.0e-30)
    route_spread = max(
        abs(left - right) for left, right in zip(final["moment"], final["virtual"])
    ) / peak
    refinement_change = _relative_l2(final["moment"], prior["moment"])
    reference_nrms = _relative_l2(reference_torque, final["moment"])
    reference_peak = max(max(abs(value) for value in reference_torque), 1.0e-30)
    peak_difference = abs(peak - reference_peak) / reference_peak
    endpoint_ratio = max(abs(final["moment"][0]), abs(final["moment"][-1])) / peak
    identity = str(summary.get("physics_identity_sha256", ""))
    reference_identity = str(reference.get("physics_identity_sha256", ""))
    monotone = all(
        left["hex_element_count"] < right["hex_element_count"]
        and left["ndof"] < right["ndof"]
        for left, right in zip(parsed_levels, parsed_levels[1:])
    )
    checks = {
        "schema_is_hdiv_hex_motor_torque_evidence": (
            summary.get("schema") == "radia.hdiv-hex-motor-torque-evidence.v1"
        ),
        "hdiv_mmm_lane_declared": summary.get("lane") == "hdiv_mmm",
        "physics_identity_is_content_addressed": (
            len(identity) == 64 and all(character in "0123456789abcdef" for character in identity)
        ),
        "all_levels_use_bdm1_project_lane_on_hex": all(
            level["space"] == _SPACE for level in parsed_levels
        ),
        "hex_elements_and_dofs_strictly_increase": monotone,
        "operator_built_once_per_level": all(
            level["operator_build_count"] == 1 for level in parsed_levels
        ),
        "moment_and_virtual_work_torque_agree": route_spread <= route_tolerance,
        "final_refinement_change_is_small": refinement_change <= refinement_tolerance,
        "torque_endpoints_are_small": endpoint_ratio <= endpoint_tolerance,
        "reference_is_executed_same_identity": (
            reference.get("executed") is True
            and identity == reference_identity
        ),
        "reference_waveform_agrees": reference_nrms <= reference_tolerance,
        "reference_peak_agrees": peak_difference <= reference_tolerance,
    }
    issues = [name for name, ok in checks.items() if not ok]
    return {
        "policy": "hdiv_hex_motor_torque_gate_v1",
        "status": "ok" if not issues else "needs_attention",
        "checks": checks,
        "issues": issues,
        "metrics": {
            "level_count": len(parsed_levels),
            "final_hex_element_count": final["hex_element_count"],
            "final_ndof": final["ndof"],
            "peak_torque_Nm": peak,
            "moment_virtual_work_spread_relative": route_spread,
            "final_refinement_relative_change": refinement_change,
            "reference_waveform_normalized_rms": reference_nrms,
            "reference_peak_relative_difference": peak_difference,
            "endpoint_to_peak_ratio": endpoint_ratio,
        },
        "discrete_space": dict(_SPACE),
        "lesson": (
            "HDiv-MMM motor torque requires the BDM1 project lane on HEX, "
            "monotone refinement, independent moment/coenergy torque, and an "
            "executed same-identity reference."
        ),
    }


__all__ = ["hdiv_hex_motor_torque_gate"]
