"""Solver-neutral gates for symmetric field-profile cross-validation."""

from __future__ import annotations

import math
import hashlib
import json
from typing import Any


def nonlinear_magnetic_refinement_energy_gate(
    summary: dict[str, Any],
    *,
    max_refinement_growth_factor: float = 1.05,
    max_finest_pair_relative_change: float = 0.05,
    max_legendre_relative_residual: float = 1.0e-8,
    min_refinement_levels: int = 3,
) -> dict[str, Any]:
    """Gate nonlinear volume observables across a material-matched h ladder.

    Every level must use a response/material order pair of ``p``/``p-1``, a
    finite positive physical permeability field, and matching field/energy
    identities. Every level must descend from one physical refinement parent
    and satisfy ``W + W* = integral(H dot B)``. Successive average-field, RMS,
    energy, and coenergy changes must contract toward the finest mesh; a single
    close result is intentionally insufficient.
    """

    if not isinstance(summary, dict):
        raise ValueError("summary must be a mapping")
    growth_limit = float(max_refinement_growth_factor)
    finest_limit = float(max_finest_pair_relative_change)
    legendre_limit = float(max_legendre_relative_residual)
    level_limit = int(min_refinement_levels)
    if not math.isfinite(growth_limit) or growth_limit < 1.0:
        raise ValueError("max_refinement_growth_factor must be finite and >= 1")
    if not math.isfinite(finest_limit) or finest_limit < 0.0:
        raise ValueError("max_finest_pair_relative_change must be finite and nonnegative")
    if not math.isfinite(legendre_limit) or legendre_limit < 0.0:
        raise ValueError("max_legendre_relative_residual must be finite and nonnegative")
    if level_limit < 3:
        raise ValueError("min_refinement_levels must be at least 3")

    raw_levels = summary.get("levels")
    if not isinstance(raw_levels, list):
        raise ValueError("levels must be a list")

    identity_keys = (
        "material_domain",
        "coordinate_system",
        "unit_system",
        "nonlinear_state_id",
        "mesh_topology_geometry_sha256",
        "bh_table_sha256",
        "material_state_identity_sha256",
        "solution_sha256",
        "refinement_parent_identity_sha256",
    )
    levels: list[dict[str, Any]] = []
    for index, raw in enumerate(raw_levels):
        if not isinstance(raw, dict):
            raise ValueError(f"levels[{index}] must be a mapping")

        def number(name: str) -> float | None:
            try:
                value = float(raw.get(name))
            except (TypeError, ValueError):
                return None
            return value if math.isfinite(value) else None

        average_raw = raw.get("average_field_T")
        try:
            average = [float(value) for value in average_raw]
        except (TypeError, ValueError):
            average = []
        if len(average) != 3 or not all(math.isfinite(value) for value in average):
            average = []
        try:
            response_order = int(raw.get("response_order"))
            material_order = int(raw.get("material_update_order"))
        except (TypeError, ValueError):
            response_order = material_order = -1
        bounds = raw.get("physical_relative_permeability_bounds")
        try:
            permeability_bounds = [float(value) for value in bounds]
        except (TypeError, ValueError):
            permeability_bounds = []
        field_identity = raw.get("field_identity")
        energy_identity = raw.get("energy_identity")
        field_identity = field_identity if isinstance(field_identity, dict) else {}
        energy_identity = energy_identity if isinstance(energy_identity, dict) else {}
        identities_match = (
            all(
                bool(str(field_identity.get(name) or "").strip())
                and str(field_identity.get(name)).strip()
                == str(energy_identity.get(name) or "").strip()
                for name in identity_keys
            )
            and str(field_identity.get("unit") or "").strip() == "T"
            and str(energy_identity.get("unit") or "").strip() == "J"
        )
        levels.append(
            {
                "mesh_size_m": number("mesh_size_m"),
                "response_order": response_order,
                "material_update_order": material_order,
                "solver_converged": raw.get("solver_converged") is True,
                "volume_m3": number("volume_m3"),
                "average_field_T": average,
                "rms_magnitude_T": number("rms_magnitude_T"),
                "magnetic_energy_J": number("magnetic_energy_J"),
                "magnetic_coenergy_J": number("magnetic_coenergy_J"),
                "h_dot_b_integral_J": number("h_dot_b_integral_J"),
                "legendre_residual_relative": number("legendre_residual_relative"),
                "physical_relative_permeability_bounds": permeability_bounds,
                "field_energy_identity_matches": identities_match,
                "refinement_parent_identity_sha256": str(
                    energy_identity.get("refinement_parent_identity_sha256") or ""
                ).strip().lower(),
            }
        )

    mesh_sizes = [level["mesh_size_m"] for level in levels]
    mesh_ordered = all(
        mesh_sizes[index] is not None
        and mesh_sizes[index + 1] is not None
        and mesh_sizes[index] > mesh_sizes[index + 1] > 0.0
        for index in range(max(0, len(mesh_sizes) - 1))
    )
    level_validity = []
    for level in levels:
        average = level["average_field_T"]
        average_magnitude = (
            math.sqrt(sum(value * value for value in average)) if average else None
        )
        rms = level["rms_magnitude_T"]
        bounds = level["physical_relative_permeability_bounds"]
        level_validity.append(
            {
                "solver_converged": level["solver_converged"],
                "response_material_orders_compatible": (
                    level["response_order"] >= 1
                    and level["material_update_order"] == level["response_order"] - 1
                ),
                "volume_positive": level["volume_m3"] is not None and level["volume_m3"] > 0.0,
                "field_finite_nonzero": average_magnitude is not None and average_magnitude > 0.0,
                "rms_consistent_with_average": (
                    rms is not None
                    and average_magnitude is not None
                    and rms + 1.0e-12 * max(rms, average_magnitude, 1.0)
                    >= average_magnitude
                ),
                "energy_finite_nonnegative": (
                    level["magnetic_energy_J"] is not None
                    and level["magnetic_energy_J"] >= 0.0
                    and level["magnetic_coenergy_J"] is not None
                    and level["magnetic_coenergy_J"] >= 0.0
                ),
                "legendre_energy_identity_satisfied": (
                    level["h_dot_b_integral_J"] is not None
                    and level["h_dot_b_integral_J"] >= 0.0
                    and level["legendre_residual_relative"] is not None
                    and level["legendre_residual_relative"] <= legendre_limit
                ),
                "physical_permeability_positive": (
                    len(bounds) == 2
                    and all(math.isfinite(value) for value in bounds)
                    and 0.0 < bounds[0] <= bounds[1]
                ),
                "field_energy_identity_matches": level["field_energy_identity_matches"],
            }
        )

    parent_digests = [level["refinement_parent_identity_sha256"] for level in levels]
    refinement_parent_consistent = bool(parent_digests) and all(
        len(value) == 64
        and all(character in "0123456789abcdef" for character in value)
        and value == parent_digests[0]
        for value in parent_digests
    )

    pair_changes: list[dict[str, float]] = []
    for coarse, fine in zip(levels, levels[1:]):
        coarse_average = coarse["average_field_T"]
        fine_average = fine["average_field_T"]
        coarse_rms = coarse["rms_magnitude_T"]
        fine_rms = fine["rms_magnitude_T"]
        coarse_energy = coarse["magnetic_energy_J"]
        fine_energy = fine["magnetic_energy_J"]
        coarse_coenergy = coarse["magnetic_coenergy_J"]
        fine_coenergy = fine["magnetic_coenergy_J"]
        if (
            not coarse_average
            or not fine_average
            or coarse_rms is None
            or fine_rms is None
            or coarse_energy is None
            or fine_energy is None
            or coarse_coenergy is None
            or fine_coenergy is None
        ):
            pair_changes.append(
                {
                    "average": math.inf,
                    "rms": math.inf,
                    "energy": math.inf,
                    "coenergy": math.inf,
                    "combined": math.inf,
                }
            )
            continue
        average_scale = max(
            math.sqrt(sum(value * value for value in fine_average)), 1.0e-300
        )
        average_change = (
            math.sqrt(
                sum(
                    (fine_value - coarse_value) ** 2
                    for coarse_value, fine_value in zip(coarse_average, fine_average)
                )
            )
            / average_scale
        )
        rms_change = abs(fine_rms - coarse_rms) / max(abs(fine_rms), 1.0e-300)
        energy_change = abs(fine_energy - coarse_energy) / max(
            abs(fine_energy), 1.0e-300
        )
        coenergy_change = abs(fine_coenergy - coarse_coenergy) / max(
            abs(fine_coenergy), 1.0e-300
        )
        pair_changes.append(
            {
                "average": average_change,
                "rms": rms_change,
                "energy": energy_change,
                "coenergy": coenergy_change,
                "combined": max(
                    average_change, rms_change, energy_change, coenergy_change
                ),
            }
        )

    changes_contract = all(
        pair_changes[index + 1]["combined"]
        <= growth_limit * pair_changes[index]["combined"]
        for index in range(max(0, len(pair_changes) - 1))
    )
    finest_change = pair_changes[-1]["combined"] if pair_changes else math.inf
    checks = {
        "refinement_levels_sufficient": len(levels) >= level_limit,
        "mesh_sizes_strictly_decrease": mesh_ordered,
        "refinement_parent_identity_consistent": refinement_parent_consistent,
        "all_levels_valid": bool(level_validity)
        and all(all(checks.values()) for checks in level_validity),
        "field_rms_energy_changes_contract": changes_contract,
        "finest_pair_is_stable": finest_change <= finest_limit,
    }
    return {
        "policy": "nonlinear_magnetic_refinement_energy_gate_v3",
        "status": "ok" if all(checks.values()) else "needs_attention",
        "checks": checks,
        "issues": [name for name, accepted in checks.items() if not accepted],
        "level_checks": level_validity,
        "pair_relative_changes": pair_changes,
        "metrics": {
            "refinement_level_count": len(levels),
            "finest_pair_relative_change": finest_change,
        },
        "tolerances": {
            "max_refinement_growth_factor": growth_limit,
            "max_finest_pair_relative_change": finest_limit,
            "min_refinement_levels": level_limit,
            "max_legendre_relative_residual": legendre_limit,
        },
        "notes": [
            "one mesh or one point cannot establish nonlinear spatial convergence",
            "field and energy observables must bind the same mesh, solution, B-H table, material state, domain, and frame",
            "every level must descend from one explicit physical refinement parent",
            "energy and coenergy must close independently against the same-state H dot B integral",
            "the gate diagnoses evidence quality and does not assert cross-solver parity by itself",
        ],
    }


def nonlinear_magnetic_field_energy_parity_gate(
    summary: dict[str, Any],
    *,
    max_average_vector_relative_difference: float = 0.05,
    max_rms_magnitude_relative_difference: float = 0.05,
    max_energy_relative_difference: float = 0.05,
    max_coenergy_relative_difference: float = 0.05,
) -> dict[str, Any]:
    """Compare nonlinear field and energy only after physical identity matches."""

    if not isinstance(summary, dict):
        raise ValueError("summary must be a mapping")
    tolerances = {
        "average_vector": float(max_average_vector_relative_difference),
        "rms_magnitude": float(max_rms_magnitude_relative_difference),
        "energy": float(max_energy_relative_difference),
        "coenergy": float(max_coenergy_relative_difference),
    }


    if any(not math.isfinite(value) or value < 0.0 for value in tolerances.values()):
        raise ValueError("all parity tolerances must be finite and nonnegative")

    identity_keys = (
        "geometry_identity_sha256",
        "material_identity_sha256",
        "excitation_identity_sha256",
        "bh_table_sha256",
        "constitutive_interpolation",
        "constitutive_extrapolation",
        "magnetic_anisotropy",
        "region_labels",
        "coordinate_system",
        "unit_system",
        "analysis_kind",
        "case_index",
        "time_semantics",
    )

    def normalized_lane(name: str) -> dict[str, Any]:
        raw = summary.get(name)
        if not isinstance(raw, dict):
            raise ValueError(f"{name} must be a mapping")
        identity = raw.get("identity")
        if not isinstance(identity, dict):
            identity = {}
        try:
            average = [float(value) for value in raw.get("average_field_T")]
        except (TypeError, ValueError):
            average = []
        if len(average) != 3 or not all(math.isfinite(value) for value in average):
            average = []

        def number(key: str) -> float | None:
            try:
                value = float(raw.get(key))
            except (TypeError, ValueError):
                return None
            return value if math.isfinite(value) else None

        return {
            "identity": identity,
            "average_field_T": average,
            "rms_magnitude_T": number("rms_magnitude_T"),
            "magnetic_energy_J": number("magnetic_energy_J"),
            "magnetic_coenergy_J": number("magnetic_coenergy_J"),
        }

    candidate = normalized_lane("candidate")
    reference = normalized_lane("reference")
    identity_checks = {
        key: bool(candidate["identity"].get(key) not in (None, "", []))
        and candidate["identity"].get(key) == reference["identity"].get(key)
        for key in identity_keys
    }
    identities_match = all(identity_checks.values())
    candidate_values_valid = (
        bool(candidate["average_field_T"])
        and candidate["rms_magnitude_T"] is not None
        and candidate["rms_magnitude_T"] >= 0.0
        and candidate["magnetic_energy_J"] is not None
        and candidate["magnetic_energy_J"] >= 0.0
        and candidate["magnetic_coenergy_J"] is not None
        and candidate["magnetic_coenergy_J"] >= 0.0
    )
    reference_values_valid = (
        bool(reference["average_field_T"])
        and reference["rms_magnitude_T"] is not None
        and reference["rms_magnitude_T"] >= 0.0
        and reference["magnetic_energy_J"] is not None
        and reference["magnetic_energy_J"] >= 0.0
        and reference["magnetic_coenergy_J"] is not None
        and reference["magnetic_coenergy_J"] >= 0.0
    )

    relative_differences: dict[str, float | None] = {
        "average_vector": None,
        "rms_magnitude": None,
        "energy": None,
        "coenergy": None,
    }
    if identities_match and candidate_values_valid and reference_values_valid:
        reference_average_norm = math.sqrt(
            sum(value * value for value in reference["average_field_T"])
        )
        relative_differences = {
            "average_vector": math.sqrt(
                sum(
                    (candidate_value - reference_value) ** 2
                    for candidate_value, reference_value in zip(
                        candidate["average_field_T"], reference["average_field_T"]
                    )
                )
            )
            / max(reference_average_norm, 1.0e-300),
            "rms_magnitude": abs(
                candidate["rms_magnitude_T"] - reference["rms_magnitude_T"]
            )
            / max(abs(reference["rms_magnitude_T"]), 1.0e-300),
            "energy": abs(
                candidate["magnetic_energy_J"] - reference["magnetic_energy_J"]
            )
            / max(abs(reference["magnetic_energy_J"]), 1.0e-300),
            "coenergy": abs(
                candidate["magnetic_coenergy_J"]
                - reference["magnetic_coenergy_J"]
            )
            / max(abs(reference["magnetic_coenergy_J"]), 1.0e-300),
        }

    checks = {
        "physical_identity_matches": identities_match,
        "candidate_observables_valid": candidate_values_valid,
        "reference_observables_valid": reference_values_valid,
        "average_vector_matches": relative_differences["average_vector"] is not None
        and relative_differences["average_vector"] <= tolerances["average_vector"],
        "rms_magnitude_matches": relative_differences["rms_magnitude"] is not None
        and relative_differences["rms_magnitude"] <= tolerances["rms_magnitude"],
        "energy_matches": relative_differences["energy"] is not None
        and relative_differences["energy"] <= tolerances["energy"],
        "coenergy_matches": relative_differences["coenergy"] is not None
        and relative_differences["coenergy"] <= tolerances["coenergy"],
    }
    return {
        "policy": "nonlinear_magnetic_field_energy_parity_gate_v1",
        "status": "ok" if all(checks.values()) else "needs_attention",
        "checks": checks,
        "issues": [name for name, accepted in checks.items() if not accepted],
        "identity_checks": identity_checks,
        "comparison_performed": identities_match
        and candidate_values_valid
        and reference_values_valid,
        "relative_differences": relative_differences,
        "tolerances": tolerances,
        "notes": [
            "numeric parity is never evaluated before geometry, material, B-H interpolation/extrapolation, excitation, region, frame, units, case, and time semantics match",
            "an accepted contract is evidence for these observables and tolerances, not universal solver equivalence",
        ],
    }


def nonlinear_constitutive_response_parity_gate(
    summary: dict[str, Any],
    *,
    max_response_relative_difference: float = 1.0e-6,
    max_integrability_relative_residual: float = 1.0e-8,
    min_response_samples: int = 5,
) -> dict[str, Any]:
    """Compare realized nonlinear material responses on one explicit H grid."""

    if not isinstance(summary, dict):
        raise ValueError("summary must be a mapping")
    response_limit = float(max_response_relative_difference)
    identity_limit = float(max_integrability_relative_residual)
    sample_limit = int(min_response_samples)
    if not math.isfinite(response_limit) or response_limit < 0.0:
        raise ValueError("max_response_relative_difference must be finite and nonnegative")
    if not math.isfinite(identity_limit) or identity_limit < 0.0:
        raise ValueError("max_integrability_relative_residual must be finite and nonnegative")
    if sample_limit < 5:
        raise ValueError("min_response_samples must be at least 5")

    identity_keys = (
        "bh_table_sha256",
        "material_model",
        "magnetic_anisotropy",
        "H_unit",
        "B_unit",
        "differential_permeability_unit",
        "energy_density_unit",
    )
    array_keys = (
        "H_A_per_m",
        "B_T",
        "differential_permeability_H_per_m",
        "energy_density_J_per_m3",
        "coenergy_density_J_per_m3",
        "d_energy_d_B_A_per_m",
        "d_coenergy_d_H_T",
    )

    def grid_digest(values: list[float]) -> str:
        encoded = json.dumps(
            values, ensure_ascii=True, separators=(",", ":")
        ).encode("ascii")
        return hashlib.sha256(encoded).hexdigest()

    def relative_residual(actual: list[float], expected: list[float]) -> float:
        numerator = math.sqrt(sum((a - b) ** 2 for a, b in zip(actual, expected)))
        denominator = max(
            math.sqrt(sum(value * value for value in expected)), 1.0e-300
        )
        return numerator / denominator

    def normalized_lane(name: str) -> dict[str, Any]:
        raw = summary.get(name)
        if not isinstance(raw, dict):
            raise ValueError(f"{name} must be a mapping")
        identity = raw.get("identity")
        if not isinstance(identity, dict):
            identity = {}
        arrays: dict[str, list[float]] = {}
        for key in array_keys:
            value = raw.get(key)
            if not isinstance(value, (list, tuple)):
                arrays[key] = []
                continue
            try:
                arrays[key] = [float(item) for item in value]
            except (TypeError, ValueError):
                arrays[key] = []
        lengths = {len(value) for value in arrays.values()}
        finite = bool(arrays["H_A_per_m"]) and all(
            math.isfinite(item) for values in arrays.values() for item in values
        )
        aligned = len(lengths) == 1
        h_values = arrays["H_A_per_m"]
        b_values = arrays["B_T"]
        tangent = arrays["differential_permeability_H_per_m"]
        monotone_h = bool(h_values) and h_values[0] >= 0.0 and all(
            right > left for left, right in zip(h_values, h_values[1:])
        )
        monotone_b = bool(b_values) and all(
            right >= left for left, right in zip(b_values, b_values[1:])
        )
        tangent_nonnegative = bool(tangent) and all(value >= 0.0 for value in tangent)
        computed_grid_digest = grid_digest(h_values) if h_values else ""
        bh_digest = str(identity.get("bh_table_sha256") or "")
        response_digest = str(identity.get("response_grid_sha256") or "")
        valid_hex = set("0123456789abcdefABCDEF")
        bh_digest_valid = len(bh_digest) == 64 and set(bh_digest) <= valid_hex
        response_digest_valid = (
            len(response_digest) == 64 and set(response_digest) <= valid_hex
        )
        residuals = {
            "legendre": math.inf,
            "energy_derivative": math.inf,
            "coenergy_derivative": math.inf,
        }
        interval_integrability = {
            "coenergy_in_monotone_bounds": False,
            "energy_in_monotone_bounds": False,
        }
        if finite and aligned and len(h_values) >= sample_limit:
            energy = arrays["energy_density_J_per_m3"]
            coenergy = arrays["coenergy_density_J_per_m3"]
            residuals = {
                "legendre": relative_residual(
                    [w + ws for w, ws in zip(energy, coenergy)],
                    [h * b for h, b in zip(h_values, b_values)],
                ),
                "energy_derivative": relative_residual(
                    arrays["d_energy_d_B_A_per_m"], h_values
                ),
                "coenergy_derivative": relative_residual(
                    arrays["d_coenergy_d_H_T"], b_values
                ),
            }
            tolerance = 1.0e-12
            coenergy_bounds = []
            energy_bounds = []
            for index in range(len(h_values) - 1):
                delta_h = h_values[index + 1] - h_values[index]
                delta_b = b_values[index + 1] - b_values[index]
                delta_coenergy = coenergy[index + 1] - coenergy[index]
                delta_energy = energy[index + 1] - energy[index]
                coenergy_bounds.append(
                    delta_h * b_values[index] - tolerance
                    <= delta_coenergy
                    <= delta_h * b_values[index + 1] + tolerance
                )
                energy_bounds.append(
                    h_values[index] * delta_b - tolerance
                    <= delta_energy
                    <= h_values[index + 1] * delta_b + tolerance
                )
            interval_integrability = {
                "coenergy_in_monotone_bounds": all(coenergy_bounds),
                "energy_in_monotone_bounds": all(energy_bounds),
            }
        checks = {
            "sample_count_sufficient": len(h_values) >= sample_limit,
            "arrays_finite_and_aligned": finite and aligned,
            "H_grid_strictly_increasing": monotone_h,
            "B_response_monotone": monotone_b,
            "differential_permeability_nonnegative": tangent_nonnegative,
            "SI_units_explicit": identity.get("H_unit") == "A/m"
            and identity.get("B_unit") == "T"
            and identity.get("differential_permeability_unit") == "H/m"
            and identity.get("energy_density_unit") == "J/m^3",
            "identity_digests_valid": bh_digest_valid and response_digest_valid,
            "reported_grid_digest_matches": bool(computed_grid_digest)
            and identity.get("response_grid_sha256") == computed_grid_digest,
            "legendre_identity_satisfied": residuals["legendre"] <= identity_limit,
            "energy_derivative_identity_satisfied": (
                residuals["energy_derivative"] <= identity_limit
            ),
            "coenergy_derivative_identity_satisfied": (
                residuals["coenergy_derivative"] <= identity_limit
            ),
            **interval_integrability,
        }
        return {
            "identity": identity,
            "arrays": arrays,
            "checks": checks,
            "integrability_relative_residuals": residuals,
            "valid": all(checks.values()),
        }

    candidate = normalized_lane("candidate")
    reference = normalized_lane("reference")
    identity_checks = {
        key: bool(candidate["identity"].get(key) not in (None, "", []))
        and candidate["identity"].get(key) == reference["identity"].get(key)
        for key in identity_keys
    }
    grids_match = (
        candidate["arrays"]["H_A_per_m"] == reference["arrays"]["H_A_per_m"]
        and candidate["identity"].get("response_grid_sha256")
        == reference["identity"].get("response_grid_sha256")
    )
    comparison_performed = (
        candidate["valid"]
        and reference["valid"]
        and all(identity_checks.values())
        and grids_match
    )
    compared_keys = (
        "B_T",
        "differential_permeability_H_per_m",
        "energy_density_J_per_m3",
        "coenergy_density_J_per_m3",
    )
    differences: dict[str, float | None] = {key: None for key in compared_keys}
    if comparison_performed:
        differences = {
            key: relative_residual(candidate["arrays"][key], reference["arrays"][key])
            for key in compared_keys
        }
    checks = {
        "candidate_response_valid": candidate["valid"],
        "reference_response_valid": reference["valid"],
        "constitutive_identity_matches": all(identity_checks.values()),
        "response_grid_matches": grids_match,
        **{
            f"{key}_matches": value is not None and value <= response_limit
            for key, value in differences.items()
        },
    }
    return {
        "policy": "nonlinear_constitutive_response_parity_gate_v1",
        "status": "ok" if all(checks.values()) else "needs_attention",
        "checks": checks,
        "issues": [name for name, accepted in checks.items() if not accepted],
        "identity_checks": identity_checks,
        "lane_checks": {
            "candidate": candidate["checks"],
            "reference": reference["checks"],
        },
        "integrability_relative_residuals": {
            "candidate": candidate["integrability_relative_residuals"],
            "reference": reference["integrability_relative_residuals"],
        },
        "comparison_performed": comparison_performed,
        "relative_differences": differences,
        "tolerances": {
            "response_relative_difference": response_limit,
            "integrability_relative_residual": identity_limit,
            "minimum_response_samples": sample_limit,
        },
        "notes": [
            "equal B-H knot digests do not establish equal realized interpolation, tangent permeability, or saturation-tail response",
            "numeric comparison requires one explicit SI H grid and independently integrable energy/coenergy responses",
            "an accepted response contract is constitutive evidence only, not whole-solver equivalence",
        ],
    }


def nonlinear_magnetic_spatial_evidence_gate(
    summary: dict[str, Any],
    *,
    max_average_vector_relative_difference: float = 0.07,
    max_rms_magnitude_relative_difference: float = 0.10,
    min_tensor_gauss_samples: int = 27,
) -> dict[str, Any]:
    """Gate nonlinear magnetic evidence using a matched volume observable.

    A converged fixed-point residual or one field sample is not sufficient. The
    response and material-update polynomial orders must be compatible, and the
    candidate must be compared with an independent result carrying the same
    observable identity, field unit, and coordinate system.
    """

    if not isinstance(summary, dict):
        raise ValueError("summary must be a mapping")
    average_limit = float(max_average_vector_relative_difference)
    rms_limit = float(max_rms_magnitude_relative_difference)
    min_samples = int(min_tensor_gauss_samples)
    if any(
        not math.isfinite(value) or value < 0.0
        for value in (average_limit, rms_limit)
    ):
        raise ValueError("relative tolerances must be finite and nonnegative")
    if min_samples < 1:
        raise ValueError("min_tensor_gauss_samples must be positive")

    def finite_number(payload: dict[str, Any], name: str) -> float | None:
        try:
            value = float(payload.get(name))
        except (TypeError, ValueError):
            return None
        return value if math.isfinite(value) else None

    def finite_vector(payload: dict[str, Any], name: str) -> list[float] | None:
        raw = payload.get(name)
        if not isinstance(raw, (list, tuple)) or len(raw) != 3:
            return None
        try:
            values = [float(value) for value in raw]
        except (TypeError, ValueError):
            return None
        return values if all(math.isfinite(value) for value in values) else None

    try:
        response_order = int(summary.get("response_order"))
        material_order = int(summary.get("material_update_order"))
        sample_count = int(summary.get("sample_count", 0))
        integration_order = int(summary.get("integration_order", 0))
    except (TypeError, ValueError):
        response_order = material_order = sample_count = integration_order = -1
    observable = str(summary.get("spatial_observable") or "").strip()
    average = finite_vector(summary, "average_field_T")
    volume = finite_number(summary, "volume_m3")
    rms = finite_number(summary, "rms_magnitude_T")
    average_magnitude = (
        math.sqrt(sum(value * value for value in average)) if average is not None else None
    )

    reference = summary.get("reference")
    reference = reference if isinstance(reference, dict) else {}
    reference_average = finite_vector(reference, "average_field_T")
    reference_rms = finite_number(reference, "rms_magnitude_T")
    reference_average_magnitude = (
        math.sqrt(sum(value * value for value in reference_average))
        if reference_average is not None
        else None
    )
    average_difference = (
        math.sqrt(
            sum(
                (candidate - expected) ** 2
                for candidate, expected in zip(average, reference_average)
            )
        )
        / reference_average_magnitude
        if average is not None
        and reference_average is not None
        and reference_average_magnitude is not None
        and reference_average_magnitude > 0.0
        else math.inf
    )
    rms_difference = (
        abs(rms - reference_rms) / reference_rms
        if rms is not None and reference_rms is not None and reference_rms > 0.0
        else math.inf
    )
    identity_fields = ("observable_id", "field_unit", "coordinate_system")
    identity_matches = all(
        bool(str(summary.get(name) or "").strip())
        and str(summary.get(name)).strip() == str(reference.get(name) or "").strip()
        for name in identity_fields
    )
    sampling_sufficient = (
        observable == "volume_integral"
        and integration_order >= max(2, 2 * max(response_order, 1))
    ) or (observable == "tensor_gauss" and sample_count >= min_samples)

    checks = {
        "nonlinear_constitutive_result": summary.get("nonlinear") is True,
        "solver_converged": summary.get("solver_converged") is True,
        "not_linear_only_evidence": summary.get("linear_reference_only") is not True,
        "response_material_orders_compatible": (
            response_order >= 1 and material_order >= max(response_order - 1, 0)
        ),
        "spatial_observable_supported": observable in {"volume_integral", "tensor_gauss"},
        "spatial_sampling_sufficient": sampling_sufficient,
        "volume_positive": volume is not None and volume > 0.0,
        "field_finite_nonzero": average_magnitude is not None and average_magnitude > 0.0,
        "rms_consistent_with_average": (
            rms is not None
            and average_magnitude is not None
            and rms + 1.0e-12 * max(rms, average_magnitude, 1.0) >= average_magnitude
        ),
        "reference_identity_matches": identity_matches,
        "average_vector_matches_reference": average_difference <= average_limit,
        "rms_magnitude_matches_reference": rms_difference <= rms_limit,
    }
    return {
        "policy": "nonlinear_magnetic_spatial_evidence_gate_v1",
        "status": "ok" if all(checks.values()) else "needs_attention",
        "checks": checks,
        "issues": [name for name, accepted in checks.items() if not accepted],
        "metrics": {
            "response_order": response_order,
            "material_update_order": material_order,
            "sample_count": sample_count,
            "integration_order": integration_order,
            "volume_m3": volume,
            "average_magnitude_T": average_magnitude,
            "rms_magnitude_T": rms,
            "average_vector_relative_difference": average_difference,
            "rms_magnitude_relative_difference": rms_difference,
        },
        "tolerances": {
            "max_average_vector_relative_difference": average_limit,
            "max_rms_magnitude_relative_difference": rms_limit,
            "min_tensor_gauss_samples": min_samples,
        },
        "notes": [
            "a center-point match and a converged nonlinear iteration are necessary but insufficient",
            "linear-source agreement cannot establish nonlinear constitutive parity",
            "response p-refinement is not evidence when the nonlinear material update remains lower order",
        ],
    }


def symmetric_complex_field_curve_gate(
    axis_positions: list[float],
    field_real: list[float],
    field_imag: list[float] | None = None,
    *,
    axis_unit: str = "m",
    field_unit: str = "A/m",
    log10_relative_residual: float,
    min_sample_count: int = 9,
    max_axis_symmetry_relative: float = 1.0e-9,
    max_field_symmetry_relative: float = 2.0e-3,
    max_log10_relative_residual: float = -8.0,
) -> dict[str, Any]:
    """Gate an origin-centered complex field curve by mirror symmetry.

    Unlike :func:`symmetric_axial_field_profile_gate`, this gate accepts even
    sample counts and does not require an analytic center value.  It is useful
    for result readers that sample a full line without placing a node exactly
    at the origin.
    """

    positions = [float(value) for value in axis_positions]
    real = [float(value) for value in field_real]
    imag = (
        [float(value) for value in field_imag]
        if field_imag is not None
        else [0.0] * len(real)
    )
    count = len(positions)
    if int(min_sample_count) < 5:
        raise ValueError("min_sample_count must be at least 5")
    if len(real) != count or len(imag) != count:
        raise ValueError("axis and field arrays must have equal length")
    if not str(axis_unit).strip() or not str(field_unit).strip():
        raise ValueError("axis_unit and field_unit must be non-empty")
    if not all(
        math.isfinite(value)
        for values in (positions, real, imag)
        for value in values
    ):
        raise ValueError("axis and field values must be finite")
    residual = float(log10_relative_residual)
    axis_tol = float(max_axis_symmetry_relative)
    field_tol = float(max_field_symmetry_relative)
    residual_limit = float(max_log10_relative_residual)
    if not math.isfinite(residual):
        raise ValueError("log10_relative_residual must be finite")
    if any(not math.isfinite(value) or value < 0.0 for value in (axis_tol, field_tol)):
        raise ValueError("relative tolerances must be finite and nonnegative")
    if not math.isfinite(residual_limit) or residual_limit >= 0.0:
        raise ValueError("max_log10_relative_residual must be finite and negative")

    axis_scale = max((abs(value) for value in positions), default=0.0)
    field_scale = max(
        (abs(complex(real[index], imag[index])) for index in range(count)),
        default=0.0,
    )
    pair_count = count // 2
    axis_symmetry_relative = (
        max(
            (abs(positions[index] + positions[-1 - index]) for index in range(pair_count)),
            default=0.0,
        )
        / axis_scale
        if axis_scale > 0.0
        else math.inf
    )
    field_symmetry_relative = (
        max(
            (
                abs(
                    complex(real[index], imag[index])
                    - complex(real[-1 - index], imag[-1 - index])
                )
                for index in range(pair_count)
            ),
            default=0.0,
        )
        / field_scale
        if field_scale > 0.0
        else math.inf
    )
    strictly_increasing = all(
        positions[index + 1] > positions[index]
        for index in range(max(0, count - 1))
    )
    if count % 2:
        center_bracketed = count > 0 and abs(positions[count // 2]) <= axis_tol * axis_scale
    else:
        center_bracketed = count >= 2 and positions[count // 2 - 1] < 0.0 < positions[count // 2]

    checks = {
        "sample_count_sufficient": count >= int(min_sample_count),
        "axis_strictly_increasing": strictly_increasing,
        "axis_straddles_origin": count > 1 and positions[0] < 0.0 < positions[-1],
        "origin_sampled_or_bracketed": center_bracketed,
        "axis_is_antisymmetric": axis_symmetry_relative <= axis_tol,
        "complex_field_nonzero": field_scale > 0.0,
        "complex_field_is_mirror_symmetric": field_symmetry_relative <= field_tol,
        "solver_residual_converged": residual <= residual_limit,
    }
    return {
        "policy": "symmetric_complex_field_curve_gate_v1",
        "status": "ok" if all(checks.values()) else "needs_attention",
        "checks": checks,
        "issues": [name for name, ok in checks.items() if not ok],
        "metrics": {
            "sample_count": count,
            "pair_count": pair_count,
            "axis_min": positions[0] if count else None,
            "axis_max": positions[-1] if count else None,
            "field_scale": field_scale,
            "axis_symmetry_relative": axis_symmetry_relative,
            "field_symmetry_relative": field_symmetry_relative,
            "log10_relative_residual": residual,
        },
        "units": {"axis": str(axis_unit), "field": str(field_unit)},
        "tolerances": {
            "min_sample_count": int(min_sample_count),
            "max_axis_symmetry_relative": axis_tol,
            "max_field_symmetry_relative": field_tol,
            "max_log10_relative_residual": residual_limit,
        },
        "notes": [
            "even sample counts are valid when the origin is bracketed by the central pair",
            "the complex field is compared directly, so magnitude-only phase errors cannot pass",
            "mirror symmetry is a validation identity, not an independent absolute-field reference",
        ],
    }


def symmetric_axial_field_profile_gate(
    axis_positions: list[float],
    axial_field: list[float],
    *,
    expected_center_field: float,
    transverse_field_1: list[float] | None = None,
    transverse_field_2: list[float] | None = None,
    min_sample_count: int = 5,
    max_center_relative_error: float = 1.0e-6,
    max_symmetry_relative: float = 1.0e-9,
    max_transverse_relative: float = 1.0e-9,
    max_axis_symmetry_relative: float = 1.0e-9,
    monotonic_relative_slack: float = 1.0e-12,
) -> dict[str, Any]:
    """Gate an odd, origin-centered axial field profile against an analytic value.

    The gate is solver-neutral. It checks the whole sampled profile rather than
    accepting a center-point match alone.
    """

    positions = [float(value) for value in axis_positions]
    axial = [float(value) for value in axial_field]
    transverse_1 = (
        [float(value) for value in transverse_field_1]
        if transverse_field_1 is not None
        else [0.0] * len(positions)
    )
    transverse_2 = (
        [float(value) for value in transverse_field_2]
        if transverse_field_2 is not None
        else [0.0] * len(positions)
    )
    expected = float(expected_center_field)
    tolerances = {
        "max_center_relative_error": float(max_center_relative_error),
        "max_symmetry_relative": float(max_symmetry_relative),
        "max_transverse_relative": float(max_transverse_relative),
        "max_axis_symmetry_relative": float(max_axis_symmetry_relative),
        "monotonic_relative_slack": float(monotonic_relative_slack),
    }
    if int(min_sample_count) < 5:
        raise ValueError("min_sample_count must be at least 5")
    if any(not math.isfinite(value) or value < 0.0 for value in tolerances.values()):
        raise ValueError("tolerances must be finite and nonnegative")
    if len(positions) != len(axial):
        raise ValueError("axis_positions and axial_field must have equal length")
    if len(transverse_1) != len(positions) or len(transverse_2) != len(positions):
        raise ValueError("transverse field arrays must match axis_positions")
    if not all(math.isfinite(value) for values in (positions, axial, transverse_1, transverse_2) for value in values):
        raise ValueError("profile values must be finite")
    if not math.isfinite(expected) or expected == 0.0:
        raise ValueError("expected_center_field must be finite and nonzero")

    count = len(positions)
    center_index = count // 2
    field_scale = max((abs(value) for value in axial), default=0.0)
    axis_scale = max((abs(value) for value in positions), default=0.0)
    center_value = axial[center_index] if count else math.nan
    center_relative_error = abs(center_value - expected) / abs(expected) if count else math.inf
    symmetry_relative = (
        max((abs(axial[i] - axial[-1 - i]) for i in range(count // 2)), default=0.0)
        / field_scale
        if field_scale > 0.0
        else math.inf
    )
    axis_symmetry_relative = (
        max((abs(positions[i] + positions[-1 - i]) for i in range(count // 2)), default=0.0)
        / axis_scale
        if axis_scale > 0.0
        else math.inf
    )
    transverse_relative = (
        max((abs(value) for value in transverse_1 + transverse_2), default=0.0) / field_scale
        if field_scale > 0.0
        else math.inf
    )
    monotonic_slack = field_scale * tolerances["monotonic_relative_slack"]
    increasing_to_center = all(
        axial[index + 1] + monotonic_slack >= axial[index]
        for index in range(center_index)
    )
    decreasing_from_center = all(
        axial[index + 1] <= axial[index] + monotonic_slack
        for index in range(center_index, max(center_index, count - 1))
    )
    strictly_increasing_axis = all(
        positions[index + 1] > positions[index] for index in range(max(0, count - 1))
    )

    checks = {
        "odd_sample_count_sufficient": count >= int(min_sample_count) and count % 2 == 1,
        "axis_strictly_increasing": strictly_increasing_axis,
        "axis_straddles_origin": count > 0 and positions[0] < 0.0 < positions[-1],
        "center_sample_at_origin": (
            count > 0
            and axis_scale > 0.0
            and abs(positions[center_index]) <= tolerances["max_axis_symmetry_relative"] * axis_scale
        ),
        "axis_is_antisymmetric": axis_symmetry_relative <= tolerances["max_axis_symmetry_relative"],
        "axial_field_nonzero": field_scale > 0.0,
        "center_matches_analytic_value": center_relative_error <= tolerances["max_center_relative_error"],
        "axial_field_is_symmetric": symmetry_relative <= tolerances["max_symmetry_relative"],
        "transverse_field_is_negligible": transverse_relative <= tolerances["max_transverse_relative"],
        "field_increases_toward_center": increasing_to_center,
        "field_decreases_from_center": decreasing_from_center,
    }
    return {
        "policy": "symmetric_axial_field_profile_gate_v1",
        "status": "ok" if all(checks.values()) else "needs_attention",
        "checks": checks,
        "issues": [name for name, ok in checks.items() if not ok],
        "metrics": {
            "sample_count": count,
            "center_index": center_index,
            "center_axis": positions[center_index] if count else None,
            "center_field": center_value if count else None,
            "expected_center_field": expected,
            "center_relative_error": center_relative_error,
            "axis_symmetry_relative": axis_symmetry_relative,
            "field_symmetry_relative": symmetry_relative,
            "transverse_relative": transverse_relative,
        },
        "tolerances": {"min_sample_count": int(min_sample_count), **tolerances},
        "notes": [
            "center agreement alone is insufficient; symmetry and both half-profiles are gated",
            "transverse components must be negligible for an axial symmetry-line comparison",
        ],
    }


def dual_formulation_symmetric_field_profile_gate(
    summary: dict[str, Any],
    *,
    max_profile_relative_difference: float = 0.01,
    max_center_relative_difference: float = 0.01,
    max_symmetry_relative: float = 0.01,
    min_sample_count: int = 21,
) -> dict[str, Any]:
    """Require two formulations to agree on a nonzero symmetric profile."""

    if not isinstance(summary, dict):
        raise ValueError("summary must be a mapping")
    tolerances = (
        max_profile_relative_difference,
        max_center_relative_difference,
        max_symmetry_relative,
    )
    if any(not math.isfinite(float(value)) or float(value) < 0.0 for value in tolerances):
        raise ValueError("relative tolerances must be finite and nonnegative")
    if int(min_sample_count) < 3:
        raise ValueError("min_sample_count must be at least 3")

    def finite_number(name: str) -> float | None:
        try:
            value = float(summary.get(name))
        except (TypeError, ValueError):
            return None
        return value if math.isfinite(value) else None

    sample_count = finite_number("sample_count")
    axis_min = finite_number("axis_min")
    axis_max = finite_number("axis_max")
    center_axis = finite_number("center_axis")
    center_a = finite_number("center_value_a")
    center_b = finite_number("center_value_b")
    profile_difference = finite_number("profile_relative_l2_difference")
    center_difference = finite_number("center_relative_difference")
    symmetry_a = finite_number("symmetry_relative_a")
    symmetry_b = finite_number("symmetry_relative_b")
    field_scale = max(abs(center_a or 0.0), abs(center_b or 0.0))
    axis_span = (axis_max - axis_min) if axis_min is not None and axis_max is not None else None

    checks = {
        "component_id_recorded": bool(str(summary.get("component_id") or "").strip()),
        "field_unit_recorded": bool(str(summary.get("field_unit") or "").strip()),
        "axis_unit_recorded": bool(str(summary.get("axis_unit") or "").strip()),
        "sample_count_sufficient": sample_count is not None and sample_count >= int(min_sample_count),
        "axis_straddles_center": axis_min is not None and axis_max is not None and axis_min < 0.0 < axis_max,
        "center_sample_near_origin": (
            center_axis is not None and axis_span is not None and axis_span > 0.0
            and abs(center_axis) <= max(1e-12, 1e-6 * axis_span)
        ),
        "center_field_nonzero": field_scale > 0.0,
        "profile_formulations_agree": (
            profile_difference is not None
            and profile_difference <= float(max_profile_relative_difference)
        ),
        "center_formulations_agree": (
            center_difference is not None
            and center_difference <= float(max_center_relative_difference)
        ),
        "formulation_a_is_symmetric": (
            symmetry_a is not None and symmetry_a <= float(max_symmetry_relative)
        ),
        "formulation_b_is_symmetric": (
            symmetry_b is not None and symmetry_b <= float(max_symmetry_relative)
        ),
    }
    return {
        "policy": "dual_formulation_symmetric_field_profile_gate",
        "status": "ok" if all(checks.values()) else "needs_attention",
        "checks": checks,
        "issues": [name for name, ok in checks.items() if not ok],
        "metrics": {
            "sample_count": sample_count,
            "axis_min": axis_min,
            "axis_max": axis_max,
            "center_axis": center_axis,
            "center_value_a": center_a,
            "center_value_b": center_b,
            "profile_relative_l2_difference": profile_difference,
            "center_relative_difference": center_difference,
            "symmetry_relative_a": symmetry_a,
            "symmetry_relative_b": symmetry_b,
        },
        "tolerances": {
            "max_profile_relative_difference": float(max_profile_relative_difference),
            "max_center_relative_difference": float(max_center_relative_difference),
            "max_symmetry_relative": float(max_symmetry_relative),
            "min_sample_count": int(min_sample_count),
        },
        "notes": [
            "agreement at one point is insufficient; compare the full profile and center",
            "each formulation must independently satisfy the expected reflection symmetry",
            "a zero field cannot pass a relative-agreement gate",
        ],
    }
