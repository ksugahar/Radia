"""Evidence gate for PM motor armature reaction on BDM1-lane HEX meshes."""

from __future__ import annotations

import math
import re
from collections.abc import Mapping, Sequence
from typing import Any


_SHA256 = re.compile(r"^[0-9a-f]{64}$")
_SPACE = {
    "family": "HDiv",
    "order": 1,
    "cell_family": "HEX",
    "project_lane": "BDM1",
    "strict_name": "tensor_product_hdiv_order1",
    "simplex_analogue": "BDM1",
}


def _finite(value: Any, name: str) -> float:
    parsed = float(value)
    if not math.isfinite(parsed):
        raise ValueError(f"{name} must be finite")
    return parsed


def _nonnegative(value: Any, name: str) -> float:
    parsed = _finite(value, name)
    if parsed < 0.0:
        raise ValueError(f"{name} must be nonnegative")
    return parsed


def _correlation(value: Any, name: str) -> float:
    parsed = _finite(value, name)
    if parsed < -1.0 or parsed > 1.0:
        raise ValueError(f"{name} must lie in [-1, 1]")
    return parsed


def pm_armature_reaction_hdiv_hex_gate(
    summary: Mapping[str, Any],
    *,
    incremental_tolerance: float = 3.0e-2,
    correlation_tolerance: float = 0.999,
    absolute_tolerance: float = 3.0e-2,
) -> dict[str, Any]:
    """Validate current-induced PM load-line response without hiding baseline error."""

    if not isinstance(summary, Mapping):
        raise ValueError("summary must be a mapping")
    levels = summary.get("levels")
    if not isinstance(levels, Sequence) or isinstance(levels, (str, bytes)):
        raise ValueError("levels must be a sequence")
    if len(levels) < 3:
        raise ValueError("levels must contain at least three refinement levels")

    parsed_levels = []
    for index, level in enumerate(levels):
        if not isinstance(level, Mapping):
            raise ValueError(f"levels[{index}] must be a mapping")
        coil = level.get("coil")
        coil = coil if isinstance(coil, Mapping) else {}
        quadrature = coil.get("quadrature")
        quadrature = list(quadrature) if isinstance(quadrature, Sequence) else []
        parsed_levels.append({
            "hex_count": int(level.get("hex_element_count", 0)),
            "space": dict(level.get("discrete_space") or {}),
            "incremental_error": _nonnegative(
                level.get("delta_B_normalized_rms_difference"),
                f"levels[{index}].delta_B_normalized_rms_difference",
            ),
            "correlation": _correlation(
                level.get("delta_B_waveform_correlation"),
                f"levels[{index}].delta_B_waveform_correlation",
            ),
            "absolute_zero_error": _nonnegative(
                level.get("zero_current_B_normalized_rms_difference"),
                f"levels[{index}].zero_current_B_normalized_rms_difference",
            ),
            "absolute_load_error": _nonnegative(
                level.get("loaded_B_normalized_rms_difference"),
                f"levels[{index}].loaded_B_normalized_rms_difference",
            ),
            "zero_knee_match": level.get("zero_current_knee_classification_match") is True,
            "loaded_knee_match": level.get("loaded_knee_classification_match") is True,
            "coil_representation": coil.get("representation"),
            "quadrature": quadrature,
            "filament_count": int(coil.get("filament_count", 0)),
        })

    identity = str(summary.get("physics_identity_sha256") or "").lower()
    reference = summary.get("reference")
    reference = reference if isinstance(reference, Mapping) else {}
    counts = [level["hex_count"] for level in parsed_levels]
    incremental_errors = [level["incremental_error"] for level in parsed_levels]
    zero_errors = [level["absolute_zero_error"] for level in parsed_levels]
    load_errors = [level["absolute_load_error"] for level in parsed_levels]

    integrity_checks = {
        "schema_is_pm_armature_reaction_evidence": summary.get("schema")
        == "radia.hdiv-hex-pm-armature-reaction-evidence.v1",
        "hdiv_mmm_lane_declared": summary.get("lane") == "hdiv_mmm",
        "physics_identity_is_content_addressed": bool(_SHA256.fullmatch(identity)),
        "reference_is_executed_same_identity": reference.get("executed") is True
        and str(reference.get("physics_identity_sha256") or "").lower() == identity,
        "all_levels_use_bdm1_project_lane_on_hex": all(
            level["space"] == _SPACE for level in parsed_levels
        ),
        "matching_hex_refinement_strictly_increases": all(
            left < right for left, right in zip(counts, counts[1:])
        ) and counts[0] > 0,
        "finite_section_coil_is_reconstructed_by_tensor_quadrature": all(
            level["coil_representation"] == "finite_section_gauss_filaments"
            and level["quadrature"] == [8, 8]
            and level["filament_count"] == 64
            for level in parsed_levels
        ),
        "scope_does_not_overclaim_retirement": summary.get(
            "research_lab_retirement_ready"
        ) is False and summary.get("product_or_market_retirement_ready") is False,
    }
    incremental_checks = {
        "incremental_error_decreases_monotonically": all(
            left > right for left, right in zip(incremental_errors, incremental_errors[1:])
        ),
        "final_incremental_error_within_tolerance": (
            incremental_errors[-1] <= incremental_tolerance
        ),
        "final_increment_waveform_correlation_is_high": (
            parsed_levels[-1]["correlation"] >= correlation_tolerance
        ),
        "final_two_levels_match_knee_classification": all(
            level["zero_knee_match"] and level["loaded_knee_match"]
            for level in parsed_levels[-2:]
        ),
    }
    absolute_checks = {
        "absolute_field_error_decreases_monotonically": all(
            left > right for left, right in zip(zero_errors, zero_errors[1:])
        ) and all(left > right for left, right in zip(load_errors, load_errors[1:])),
        "final_absolute_field_within_tolerance": (
            max(zero_errors[-1], load_errors[-1]) <= absolute_tolerance
        ),
    }
    integrity_ok = all(integrity_checks.values())
    incremental_ok = integrity_ok and all(incremental_checks.values())
    absolute_ok = incremental_ok and all(absolute_checks.values())
    all_checks = {**integrity_checks, **incremental_checks, **absolute_checks}
    if absolute_ok:
        status = "validated"
    elif incremental_ok:
        status = "validated_partial"
    else:
        status = "needs_attention"
    return {
        "policy": "pm_armature_reaction_hdiv_hex_gate_v1",
        "status": status,
        "checks": all_checks,
        "issues": [name for name, passed in all_checks.items() if not passed],
        "armature_reaction_increment_validated": incremental_ok,
        "absolute_self_demagnetizing_field_validated": absolute_ok,
        "irreversible_demagnetization_state_validated": False,
        "research_lab_retirement_ready": False,
        "metrics": {
            "level_count": len(parsed_levels),
            "final_hex_element_count": counts[-1],
            "final_incremental_normalized_rms_difference": incremental_errors[-1],
            "final_increment_waveform_correlation": parsed_levels[-1]["correlation"],
            "final_absolute_normalized_rms_difference": max(
                zero_errors[-1], load_errors[-1]
            ),
        },
        "discrete_space": dict(_SPACE),
        "lesson": (
            "Validate the current-induced PM load-line increment separately from the "
            "absolute self-demagnetizing baseline. Incremental agreement alone does "
            "not validate irreversible demagnetization or complete motor retirement."
        ),
    }


def pm_absolute_demag_three_way_gate(
    summary: Mapping[str, Any],
    *,
    absolute_tolerance: float = 3.0e-2,
    fem_order_stability_tolerance: float = 5.0e-3,
    fem_outer_stability_tolerance: float = 3.0e-3,
    iteration_tolerance: float = 1.0e-6,
) -> dict[str, Any]:
    """Attribute absolute PM self-demag error using HDiv orders and H1 FEM."""

    if not isinstance(summary, Mapping):
        raise ValueError("summary must be a mapping")
    identity = str(summary.get("physics_identity_sha256") or "").lower()
    reference = summary.get("reference")
    reference = reference if isinstance(reference, Mapping) else {}
    hdiv = summary.get("hdiv")
    hdiv = hdiv if isinstance(hdiv, Mapping) else {}
    fem = summary.get("independent_fem")
    fem = fem if isinstance(fem, Mapping) else {}

    bdm1_errors = [
        _nonnegative(value, f"hdiv.bdm1_source_nrms[{index}]")
        for index, value in enumerate(hdiv.get("bdm1_source_nrms") or [])
    ]
    bdm2_errors = [
        _nonnegative(value, f"hdiv.bdm2_source_nrms[{index}]")
        for index, value in enumerate(hdiv.get("bdm2_source_nrms") or [])
    ]
    fem_outer_errors = [
        _nonnegative(value, f"independent_fem.source_nrms_outer_p2[{index}]")
        for index, value in enumerate(fem.get("source_nrms_outer_p2") or [])
    ]
    if len(bdm1_errors) < 3 or len(bdm2_errors) < 3:
        raise ValueError("hdiv BDM1 and BDM2 errors require at least three levels")
    if len(fem_outer_errors) < 3:
        raise ValueError("independent FEM outer-radius errors require at least three levels")

    fem_p3_error = _nonnegative(
        fem.get("source_nrms_p3"), "independent_fem.source_nrms_p3"
    )
    fem_p2_p3 = _nonnegative(
        fem.get("p2_p3_same_mesh_nrms"),
        "independent_fem.p2_p3_same_mesh_nrms",
    )
    bdm2_vs_fem = _nonnegative(
        hdiv.get("bdm2_final_vs_fem_nrms"),
        "hdiv.bdm2_final_vs_fem_nrms",
    )
    relative_step = _nonnegative(
        hdiv.get("bdm2_final_relative_step"),
        "hdiv.bdm2_final_relative_step",
    )
    outer_spread = max(fem_outer_errors) - min(fem_outer_errors)

    integrity_checks = {
        "schema_is_three_way_absolute_demag_evidence": summary.get("schema")
        == "radia.pm-ring-absolute-demag-attribution-evidence.v1",
        "physics_identity_is_content_addressed": bool(_SHA256.fullmatch(identity)),
        "reference_is_executed_same_identity": reference.get("executed") is True
        and str(reference.get("physics_identity_sha256") or "").lower() == identity,
        "hdiv_uses_hex_bdm1_and_bdm2": hdiv.get("cell_family") == "HEX"
        and list(hdiv.get("project_lanes") or []) == ["BDM1", "BDM2"],
        "fem_is_independent_continuous_h1_tet": fem.get("formulation")
        == "continuous_H1_scalar_potential"
        and fem.get("cell_family") == "TET",
        "fem_sector_geometry_identity_checked": fem.get("equal_sector_volumes") is True,
        "scope_does_not_overclaim_retirement": summary.get(
            "research_lab_retirement_ready"
        ) is False and summary.get("product_or_market_retirement_ready") is False,
    }
    convergence_checks = {
        "bdm1_error_decreases_with_refinement": all(
            left > right for left, right in zip(bdm1_errors, bdm1_errors[1:])
        ),
        "bdm2_error_decreases_with_refinement": all(
            left > right for left, right in zip(bdm2_errors, bdm2_errors[1:])
        ),
        "bdm2_improves_final_bdm1_error": bdm2_errors[-1] < bdm1_errors[-1],
        "hdiv_block_iteration_is_converged": relative_step <= iteration_tolerance,
        "fem_outer_boundary_is_stable": outer_spread <= fem_outer_stability_tolerance,
        "fem_order_is_stable": fem_p2_p3 <= fem_order_stability_tolerance,
        "source_is_corroborated_by_independent_fem": fem_p3_error <= absolute_tolerance,
    }
    integrity_ok = all(integrity_checks.values())
    attribution_ok = integrity_ok and all(convergence_checks.values())
    absolute_ok = attribution_ok and max(bdm2_errors[-1], bdm2_vs_fem) <= absolute_tolerance
    all_checks = {**integrity_checks, **convergence_checks}
    if not attribution_ok:
        status = "needs_attention"
        attribution = "unresolved"
    elif absolute_ok:
        status = "validated"
        attribution = "three_way_agreement"
    else:
        status = "attributed_pending"
        attribution = "hdiv_resolution_or_cross_body_operator_gap"

    return {
        "policy": "pm_absolute_demag_three_way_gate_v1",
        "status": status,
        "checks": all_checks,
        "issues": [name for name, passed in all_checks.items() if not passed],
        "attribution": attribution,
        "independent_fem_corroborates_reference": attribution_ok,
        "absolute_self_demagnetizing_field_validated": absolute_ok,
        "armature_reaction_increment_validated_by_this_gate": False,
        "research_lab_retirement_ready": False,
        "metrics": {
            "final_bdm1_source_nrms": bdm1_errors[-1],
            "final_bdm2_source_nrms": bdm2_errors[-1],
            "final_bdm2_fem_nrms": bdm2_vs_fem,
            "fem_source_p3_nrms": fem_p3_error,
            "fem_p2_p3_nrms": fem_p2_p3,
            "fem_outer_radius_nrms_spread": outer_spread,
            "hdiv_final_relative_step": relative_step,
        },
        "next_gate": (
            "Compare a monolithic or product-space HDiv coupling with the current "
            "cross-body field operator, including interface-near field and surface-"
            "charge observables, before claiming absolute self-demag validation."
        ),
        "lesson": (
            "A converged block iteration is not a mesh-convergence proof. Use an "
            "independent continuous H1 FEM with outer-boundary and p-order controls "
            "to localize segmented-PM absolute-field discrepancies."
        ),
    }


__all__ = [
    "pm_absolute_demag_three_way_gate",
    "pm_armature_reaction_hdiv_hex_gate",
]
