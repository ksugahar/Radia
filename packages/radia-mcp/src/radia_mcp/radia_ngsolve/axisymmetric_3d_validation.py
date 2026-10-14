"""Axisymmetric-reference gates for validating 3-D force results.

The helpers here are deliberately solver-independent.  A trusted
axisymmetric run already includes the full ``2*pi*r`` volume/surface weight,
so its axial force is a full 3-D revolution quantity.  A 3-D solver can use
that number as a compact validation oracle: the full-revolution axial
component should match and the transverse components should cancel.
"""
from __future__ import annotations

import math
from typing import Iterable


_AXIS_TO_INDEX = {"x": 0, "y": 1, "z": 2}


def _axis_index(axis: str) -> int:
    key = str(axis).strip().lower()
    if key not in _AXIS_TO_INDEX:
        raise ValueError("axial_axis must be one of x, y, or z")
    return _AXIS_TO_INDEX[key]


def _vector3(values: Iterable[float]) -> list[float]:
    vector = [float(value) for value in values]
    if len(vector) != 3:
        raise ValueError("force_vector_N must have exactly three components")
    if not all(math.isfinite(value) for value in vector):
        raise ValueError("force_vector_N components must be finite")
    return vector


def _rel_error(computed: float, reference: float) -> float:
    return abs(computed - reference) / max(abs(computed), abs(reference), 1.0e-300)


def axisymmetric_to_3d_force_gate(
    axisymmetric_axial_force_N: float,
    force_vector_N: Iterable[float],
    *,
    case_id: str = "axisymmetric_to_3d_force",
    axial_axis: str = "z",
    result_basis: str = "full_revolution",
    sector_angle_deg: float = 360.0,
    axial_rtol: float = 0.03,
    axial_atol_N: float = 0.0,
    transverse_rtol: float = 0.02,
    transverse_atol_N: float = 0.0,
    check_transverse_cancellation: bool | None = None,
    metadata: dict | None = None,
) -> dict:
    """Check a 3-D force vector against a full-revolution axisymmetric result.

    Parameters
    ----------
    axisymmetric_axial_force_N:
        Full 3-D axial force from an axisymmetric calculation, e.g.
        ``eggshell_force_axi`` with the ``2*pi*r`` weight.
    force_vector_N:
        3-D force vector in global Cartesian components.
    result_basis:
        ``"full_revolution"`` if ``force_vector_N`` already covers the whole
        360 degree model.  ``"symmetry_sector"`` if the vector is from a
        periodic sector and should be scaled by ``360/sector_angle_deg`` for
        the axial comparison.  Transverse cancellation is only checked by
        default for the full-revolution basis.
    """

    ref = float(axisymmetric_axial_force_N)
    if not math.isfinite(ref):
        raise ValueError("axisymmetric_axial_force_N must be finite")
    axis_i = _axis_index(axial_axis)
    force = _vector3(force_vector_N)
    basis = str(result_basis).strip().lower()
    angle = float(sector_angle_deg)
    rel_tol = float(axial_rtol)
    abs_tol = float(axial_atol_N)
    trans_rel_tol = float(transverse_rtol)
    trans_abs_tol = float(transverse_atol_N)
    if angle <= 0.0 or angle > 360.0:
        raise ValueError("sector_angle_deg must be in (0, 360]")
    if rel_tol < 0.0 or abs_tol < 0.0:
        raise ValueError("axial tolerances must be non-negative")
    if trans_rel_tol < 0.0 or trans_abs_tol < 0.0:
        raise ValueError("transverse tolerances must be non-negative")

    if basis == "full_revolution":
        scale = 1.0
        if check_transverse_cancellation is None:
            check_transverse_cancellation = True
    elif basis == "symmetry_sector":
        scale = 360.0 / angle
        if check_transverse_cancellation is None:
            check_transverse_cancellation = False
    else:
        raise ValueError("result_basis must be full_revolution or symmetry_sector")

    full_force = [scale * value for value in force]
    axial = full_force[axis_i]
    transverse = [value for i, value in enumerate(full_force) if i != axis_i]
    transverse_mag = math.sqrt(sum(value * value for value in transverse))
    axial_abs_error = abs(axial - ref)
    axial_rel_error = _rel_error(axial, ref)
    axial_allowed = max(abs_tol, rel_tol * max(abs(ref), 1.0e-300))
    reference_scale = max(abs(ref), abs(axial), 1.0e-300)
    transverse_allowed = max(trans_abs_tol, trans_rel_tol * reference_scale)

    checks = {
        "axisymmetric_reference_finite": math.isfinite(ref),
        "sector_angle_valid": 0.0 < angle <= 360.0,
        "axial_component_matches_axisymmetric_reference": axial_abs_error <= axial_allowed,
    }
    if check_transverse_cancellation:
        checks["transverse_components_cancel"] = transverse_mag <= transverse_allowed
    else:
        checks["transverse_components_cancel"] = "not_checked_for_symmetry_sector"

    ok = all(value is True for value in checks.values() if value != "not_checked_for_symmetry_sector")
    result = {
        "policy": "axisymmetric_to_3d_force_gate",
        "case_id": str(case_id),
        "status": "ok" if ok else "needs_attention",
        "axisymmetric_reference": {
            "axial_axis": str(axial_axis).strip().lower(),
            "axial_force_N": ref,
            "quantity_basis": "full_3d_revolution_from_2pi_r_weight",
        },
        "three_d_result": {
            "input_force_vector_N": force,
            "result_basis": basis,
            "sector_angle_deg": angle,
            "scale_to_full_revolution": scale,
            "full_revolution_force_vector_N": full_force,
            "axial_component_N": axial,
            "transverse_components_N": transverse,
            "transverse_magnitude_N": transverse_mag,
        },
        "errors": {
            "axial_abs_error_N": axial_abs_error,
            "axial_rel_error": axial_rel_error,
            "axial_allowed_error_N": axial_allowed,
            "transverse_magnitude_N": transverse_mag,
            "transverse_allowed_magnitude_N": transverse_allowed,
        },
        "tolerances": {
            "axial_rtol": rel_tol,
            "axial_atol_N": abs_tol,
            "transverse_rtol": trans_rel_tol,
            "transverse_atol_N": trans_abs_tol,
        },
        "checks": checks,
        "required_result_contract": [
            "axisymmetric result must state that it is full 3D, not per-radian",
            "3D result must state full_revolution or symmetry_sector basis",
            "force vector must be global Cartesian components",
            "axis direction and sign convention must be recorded",
            "mesh/order/solver versions and timing should be saved with the result JSON",
        ],
    }
    if metadata:
        result["metadata"] = dict(metadata)
    return result


def axisymmetric_to_3d_validation_plan(
    case_id: str,
    *,
    axial_axis: str = "z",
    preferred_3d_route: str = "revolved_vol_or_occ_3d_model",
    validation_root: str = "validation_test/force_validation",
) -> dict:
    """Return the reusable validation plan for axisymmetric-to-3D checks."""

    axis = str(axial_axis).strip().lower()
    _axis_index(axis)
    root = str(validation_root).strip()
    return {
        "policy": "axisymmetric_to_3d_validation_plan",
        "case_id": str(case_id),
        "route": preferred_3d_route,
        "reference_quantity": {
            "name": "axisymmetric_axial_force_N",
            "basis": "full_3d_revolution_from_2pi_r_weight",
            "recommended_extractor": "radia_mcp.radia_ngsolve.force.eggshell_force_axi",
        },
        "three_d_quantity": {
            "name": "force_vector_N",
            "basis": "full_revolution preferred; symmetry_sector allowed for axial check",
            "recommended_extractor": "radia_mcp.radia_ngsolve.force.eggshell_force",
            "axial_axis": axis,
        },
        "required_artifacts": [
            f"{root}/validation_axisymmetric_to_3d_force_gate_summary.json",
            "axisymmetric mesh/result artifact with solver version and 2*pi*r convention",
            "3D mesh/result artifact with force-vector basis, sector angle, and timing",
            "material/current convention manifest shared by the axisymmetric and 3D runs",
        ],
        "required_checks": [
            "axial 3D force matches the axisymmetric full-revolution reference",
            "full-revolution transverse force components cancel",
            "3D sector outputs record sector angle before any scaling",
            "mesh refinement or shell-band sweep is recorded for the 3D force extractor",
        ],
        "recommended_tolerances": {
            "axial_rtol": 0.03,
            "transverse_rtol": 0.02,
            "mesh_refinement_target_rtol": 0.02,
        },
    }


def magnetic_material_pair_force_gate(
    reference_force_vector_N: Iterable[float],
    computed_force_vector_N: Iterable[float],
    *,
    case_id: str = "magnetic_material_pair_force",
    axial_axis: str = "z",
    reference_method: str = "dipole_dipole_or_coenergy",
    axial_rtol: float = 0.03,
    axial_atol_N: float = 0.0,
    vector_rtol: float = 0.03,
    vector_atol_N: float = 0.0,
    transverse_rtol: float = 0.005,
    transverse_atol_N: float = 0.0,
    metadata: dict | None = None,
) -> dict:
    """Check a 3-D magnetic-material force vector against a reference vector.

    This gate is for forces between magnetized or magnetizable bodies, where
    ``J x B`` is not the right observable.  The reference may come from
    coenergy/virtual work, weighted Maxwell stress, or an analytic magnetic
    dipole pair.  The gate keeps the vector basis explicit and separately
    checks the axial component and transverse cancellation/match.
    """

    axis_i = _axis_index(axial_axis)
    reference = _vector3(reference_force_vector_N)
    computed = _vector3(computed_force_vector_N)
    axial_rel_tol = float(axial_rtol)
    axial_abs_tol = float(axial_atol_N)
    vec_rel_tol = float(vector_rtol)
    vec_abs_tol = float(vector_atol_N)
    trans_rel_tol = float(transverse_rtol)
    trans_abs_tol = float(transverse_atol_N)
    if min(axial_rel_tol, axial_abs_tol, vec_rel_tol, vec_abs_tol, trans_rel_tol, trans_abs_tol) < 0.0:
        raise ValueError("force tolerances must be non-negative")

    diff = [computed[i] - reference[i] for i in range(3)]
    ref_norm = math.sqrt(sum(value * value for value in reference))
    comp_norm = math.sqrt(sum(value * value for value in computed))
    diff_norm = math.sqrt(sum(value * value for value in diff))
    scale = max(ref_norm, comp_norm, 1.0e-300)
    vector_rel_error = diff_norm / scale
    vector_allowed = max(vec_abs_tol, vec_rel_tol * max(ref_norm, 1.0e-300))

    axial_error = abs(computed[axis_i] - reference[axis_i])
    axial_rel_error = _rel_error(computed[axis_i], reference[axis_i])
    axial_allowed = max(axial_abs_tol, axial_rel_tol * max(abs(reference[axis_i]), 1.0e-300))

    transverse_reference = [value for i, value in enumerate(reference) if i != axis_i]
    transverse_computed = [value for i, value in enumerate(computed) if i != axis_i]
    transverse_diff = [
        transverse_computed[i] - transverse_reference[i]
        for i in range(len(transverse_reference))
    ]
    transverse_reference_mag = math.sqrt(sum(value * value for value in transverse_reference))
    transverse_computed_mag = math.sqrt(sum(value * value for value in transverse_computed))
    transverse_error = math.sqrt(sum(value * value for value in transverse_diff))
    transverse_allowed = max(trans_abs_tol, trans_rel_tol * scale)

    checks = {
        "reference_force_vector_finite": all(math.isfinite(value) for value in reference),
        "computed_force_vector_finite": all(math.isfinite(value) for value in computed),
        "axial_component_matches_reference": axial_error <= axial_allowed,
        "force_vector_matches_reference": diff_norm <= vector_allowed,
        "transverse_components_match_reference": transverse_error <= transverse_allowed,
    }
    ok = all(value is True for value in checks.values())
    result = {
        "policy": "magnetic_material_pair_force_gate",
        "case_id": str(case_id),
        "status": "ok" if ok else "needs_attention",
        "reference_method": str(reference_method),
        "axis": {
            "axial_axis": str(axial_axis).strip().lower(),
            "force_basis": "global_cartesian_components",
        },
        "reference": {
            "force_vector_N": reference,
            "force_norm_N": ref_norm,
            "transverse_components_N": transverse_reference,
            "transverse_magnitude_N": transverse_reference_mag,
        },
        "computed": {
            "force_vector_N": computed,
            "force_norm_N": comp_norm,
            "transverse_components_N": transverse_computed,
            "transverse_magnitude_N": transverse_computed_mag,
        },
        "errors": {
            "component_abs_error_N": [abs(value) for value in diff],
            "vector_abs_error_N": diff_norm,
            "vector_rel_error": vector_rel_error,
            "vector_allowed_error_N": vector_allowed,
            "axial_abs_error_N": axial_error,
            "axial_rel_error": axial_rel_error,
            "axial_allowed_error_N": axial_allowed,
            "transverse_abs_error_N": transverse_error,
            "transverse_allowed_error_N": transverse_allowed,
        },
        "tolerances": {
            "axial_rtol": axial_rel_tol,
            "axial_atol_N": axial_abs_tol,
            "vector_rtol": vec_rel_tol,
            "vector_atol_N": vec_abs_tol,
            "transverse_rtol": trans_rel_tol,
            "transverse_atol_N": trans_abs_tol,
        },
        "checks": checks,
        "required_result_contract": [
            "state whether the force is magnetic-material/Maxwell-stress/coenergy, not Lorentz JxB",
            "record source and target material selections",
            "record global Cartesian force-vector basis and sign convention",
            "record mesh/order/solver versions and timing with the result JSON",
        ],
    }
    if metadata:
        result["metadata"] = dict(metadata)
    return result
