"""Solver-neutral gates for symmetric field-profile cross-validation."""

from __future__ import annotations

import math
import hashlib
import json
from typing import Any


def build_constitutive_comparison_candidate(
    bh_table: list[list[float]],
    h_values: list[float],
    *,
    constitutive_interpolation: str = "monotone_pchip",
) -> dict[str, Any]:
    """Build a comparison response without changing Radia's production law."""

    mode = str(constitutive_interpolation)
    if mode not in {"monotone_pchip", "piecewise_linear"}:
        raise ValueError(
            "constitutive interpolation must be monotone_pchip or piecewise_linear"
        )
    if not isinstance(bh_table, list) or len(bh_table) < 2:
        raise ValueError("bh_table must contain at least two [H, B] rows")
    try:
        table = [[float(row[0]), float(row[1])] for row in bh_table]
        grid = [float(value) for value in h_values]
    except (IndexError, TypeError, ValueError) as exc:
        raise ValueError("bh_table and h_values must contain numeric values") from exc
    if any(not math.isfinite(value) for row in table for value in row):
        raise ValueError("bh_table must contain finite H and B values")
    h_tab = [row[0] for row in table]
    b_tab = [row[1] for row in table]
    if any(right <= left for left, right in zip(h_tab, h_tab[1:])):
        raise ValueError("bh_table H values must be strictly increasing")
    if any(right < left for left, right in zip(b_tab, b_tab[1:])):
        raise ValueError("bh_table B values must be non-decreasing")
    if len(grid) < 2:
        raise ValueError("h_values must contain at least two samples")
    if any(not math.isfinite(value) or value < 0.0 for value in grid):
        raise ValueError("h_values must be finite and nonnegative")
    if any(right <= left for left, right in zip(grid, grid[1:])):
        raise ValueError("h_values must be strictly increasing")

    mu_0 = 4.0 * math.pi * 1.0e-7
    widths = [right - left for left, right in zip(h_tab, h_tab[1:])]
    slopes = [
        (right - left) / width
        for left, right, width in zip(b_tab, b_tab[1:], widths)
    ]

    def endpoint_slope(h0, h1, d0, d1):
        value = ((2.0 * h0 + h1) * d0 - h0 * d1) / (h0 + h1)
        if value * d0 <= 0.0:
            return 0.0
        if d0 * d1 < 0.0 and abs(value) > abs(3.0 * d0):
            return 3.0 * d0
        return value

    if len(table) == 2:
        pchip_derivatives = [slopes[0], slopes[0]]
    else:
        pchip_derivatives = [
            endpoint_slope(widths[0], widths[1], slopes[0], slopes[1])
        ]
        for index in range(1, len(table) - 1):
            left_slope = slopes[index - 1]
            right_slope = slopes[index]
            if left_slope * right_slope <= 0.0:
                pchip_derivatives.append(0.0)
            else:
                left_weight = 2.0 * widths[index] + widths[index - 1]
                right_weight = widths[index] + 2.0 * widths[index - 1]
                pchip_derivatives.append(
                    (left_weight + right_weight)
                    / (left_weight / left_slope + right_weight / right_slope)
                )
        pchip_derivatives.append(
            endpoint_slope(widths[-1], widths[-2], slopes[-1], slopes[-2])
        )

    pchip_coefficients = []
    for index, secant in enumerate(slopes):
        width = widths[index]
        left_derivative = pchip_derivatives[index]
        right_derivative = pchip_derivatives[index + 1]
        pchip_coefficients.append(
            (
                float(b_tab[index]),
                left_derivative,
                (3.0 * secant - 2.0 * left_derivative - right_derivative) / width,
                (left_derivative + right_derivative - 2.0 * secant) / width**2,
            )
        )
    cumulative_linear = [float(b_tab[0] * h_tab[0])]
    for index, slope in enumerate(slopes):
        delta = float(h_tab[index + 1] - h_tab[index])
        cumulative_linear.append(
            cumulative_linear[-1]
            + float(b_tab[index]) * delta
            + 0.5 * float(slope) * delta**2
        )
    cumulative_pchip = [float(b_tab[0] * h_tab[0])]
    for width, coefficients in zip(widths, pchip_coefficients):
        c0, c1, c2, c3 = coefficients
        cumulative_pchip.append(
            cumulative_pchip[-1]
            + c0 * width
            + 0.5 * c1 * width**2
            + c2 * width**3 / 3.0
            + 0.25 * c3 * width**4
        )

    b_values = []
    tangent_values = []
    coenergy_values = []
    for raw_h_value in grid:
        h_value = float(raw_h_value)
        if h_value < h_tab[0]:
            b_value = float(b_tab[0])
            tangent = 0.0
            coenergy = float(b_tab[0]) * h_value
        elif h_value > h_tab[-1]:
            delta = h_value - float(h_tab[-1])
            b_value = float(b_tab[-1]) + mu_0 * delta
            tangent = mu_0
            value_at_max = (
                cumulative_pchip[-1]
                if mode == "monotone_pchip"
                else cumulative_linear[-1]
            )
            coenergy = (
                value_at_max
                + float(b_tab[-1]) * delta
                + 0.5 * mu_0 * delta**2
            )
        else:
            interval = next(
                (
                    index
                    for index in range(len(slopes))
                    if h_tab[index] <= h_value < h_tab[index + 1]
                ),
                len(slopes) - 1,
            )
            delta = h_value - float(h_tab[interval])
            if mode == "monotone_pchip":
                c0, c1, c2, c3 = pchip_coefficients[interval]
                b_value = c0 + c1 * delta + c2 * delta**2 + c3 * delta**3
                tangent = c1 + 2.0 * c2 * delta + 3.0 * c3 * delta**2
                coenergy = (
                    cumulative_pchip[interval]
                    + c0 * delta
                    + 0.5 * c1 * delta**2
                    + c2 * delta**3 / 3.0
                    + 0.25 * c3 * delta**4
                )
            else:
                tangent = float(slopes[interval])
                b_value = float(b_tab[interval]) + tangent * delta
                coenergy = (
                    cumulative_linear[interval]
                    + float(b_tab[interval]) * delta
                    + 0.5 * tangent * delta**2
                )
        b_values.append(float(b_value))
        tangent_values.append(float(tangent))
        coenergy_values.append(float(coenergy))

    canonical_table = table
    canonical_grid = grid

    def digest(value: object) -> str:
        encoded = json.dumps(
            value, ensure_ascii=True, separators=(",", ":")
        ).encode("ascii")
        return hashlib.sha256(encoded).hexdigest()

    return {
        "schema": "radia.nonlinear-constitutive-comparison-candidate.v1",
        "identity": {
            "bh_table_sha256": digest(canonical_table),
            "response_grid_sha256": digest(canonical_grid),
            "material_model": "single_valued_isotropic_soft_magnetic",
            "magnetic_anisotropy": "isotropic",
            "constitutive_interpolation": mode,
            "constitutive_extrapolation": "vacuum_slope",
            "candidate_only": True,
            "solver_runtime_mode": mode == "monotone_pchip",
            "radia_production_interpolation": "monotone_pchip",
            "H_unit": "A/m",
            "B_unit": "T",
            "differential_permeability_unit": "H/m",
            "energy_density_unit": "J/m^3",
        },
        "H_A_per_m": canonical_grid,
        "B_T": b_values,
        "differential_permeability_H_per_m": tangent_values,
        "energy_density_J_per_m3": [
            h_value * b_value - coenergy
            for h_value, b_value, coenergy in zip(
                canonical_grid, b_values, coenergy_values
            )
        ],
        "coenergy_density_J_per_m3": coenergy_values,
        "d_energy_d_B_A_per_m": canonical_grid,
        "d_coenergy_d_H_T": b_values,
        "last_table_H_A_per_m": h_tab[-1],
        "last_table_B_T": b_tab[-1],
        "claim_boundary": (
            "The piecewise-linear option is a comparison candidate only; Radia "
            "and NGSolve production solves retain the monotone-PCHIP law."
        ),
    }


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


def nonlinear_field_energy_identity_gate_v5(
    summary: dict[str, Any],
    *,
    max_energy_derivative_relative_residual: float = 1.0e-8,
    min_response_samples: int = 4,
) -> dict[str, Any]:
    """Validate canonical geometry/frame/refinement identity and energy derivatives.

    This is deliberately solver-neutral. It accepts no numeric parity claim until
    both lanes carry recomputable geometry and frame identities, a monotone unique
    refinement ladder, and discrete constitutive energy identities.
    """

    if not isinstance(summary, dict):
        raise ValueError("summary must be a mapping")
    residual_limit = float(max_energy_derivative_relative_residual)
    sample_limit = int(min_response_samples)
    if not math.isfinite(residual_limit) or residual_limit < 0.0:
        raise ValueError(
            "max_energy_derivative_relative_residual must be finite and nonnegative"
        )
    if sample_limit < 4:
        raise ValueError("min_response_samples must be at least 4")

    def digest(value: object) -> str:
        encoded = json.dumps(
            value, ensure_ascii=True, sort_keys=True, separators=(",", ":")
        ).encode("utf-8")
        return hashlib.sha256(encoded).hexdigest()

    def is_sha256(value: object) -> bool:
        text = str(value or "").strip().lower()
        return len(text) == 64 and all(
            character in "0123456789abcdef" for character in text
        )

    def frame_ok(matrix: object) -> bool:
        if not isinstance(matrix, list) or len(matrix) != 3:
            return False
        if any(
            not isinstance(row, list)
            or len(row) != 3
            or any(not isinstance(value, (int, float)) or not math.isfinite(float(value)) for value in row)
            for row in matrix
        ):
            return False
        rows = [[float(value) for value in row] for row in matrix]
        orthogonal = all(
            abs(
                sum(rows[i][axis] * rows[j][axis] for axis in range(3))
                - (1.0 if i == j else 0.0)
            )
            <= 1.0e-10
            for i in range(3)
            for j in range(3)
        )
        determinant = (
            rows[0][0] * (rows[1][1] * rows[2][2] - rows[1][2] * rows[2][1])
            - rows[0][1] * (rows[1][0] * rows[2][2] - rows[1][2] * rows[2][0])
            + rows[0][2] * (rows[1][0] * rows[2][1] - rows[1][1] * rows[2][0])
        )
        return orthogonal and determinant > 0.0

    def refinement_signature(identity: dict[str, Any]) -> tuple[bool, list[object]]:
        raw_levels = identity.get("refinement_levels")
        if not isinstance(raw_levels, list) or not raw_levels:
            return False, []
        level_ids: list[str] = []
        counts: list[int] = []
        mesh_digests: list[str] = []
        for level in raw_levels:
            if not isinstance(level, dict):
                return False, []
            level_id = str(level.get("level_id") or "").strip()
            digest_value = str(level.get("mesh_identity_sha256") or "").strip().lower()
            try:
                element_count = int(level.get("element_count"))
            except (TypeError, ValueError):
                return False, []
            if not level_id or element_count <= 0 or not is_sha256(digest_value):
                return False, []
            level_ids.append(level_id)
            counts.append(element_count)
            mesh_digests.append(digest_value)
        valid = (
            len(level_ids) == len(set(level_ids))
            and len(mesh_digests) == len(set(mesh_digests))
            and all(left < right for left, right in zip(counts, counts[1:]))
        )
        return valid, [level_ids, counts, mesh_digests]

    def lane_identity(name: str) -> dict[str, Any]:
        raw = summary.get(name)
        if not isinstance(raw, dict):
            raise ValueError(f"{name} must be a mapping")
        identity = raw.get("identity")
        if not isinstance(identity, dict):
            identity = {}
        geometry = identity.get("canonical_geometry")
        frame = identity.get("coordinate_frame_matrix")
        geometry_digest = str(identity.get("geometry_canonical_sha256") or "").strip().lower()
        frame_digest = str(identity.get("coordinate_frame_sha256") or "").strip().lower()
        response = raw.get("constitutive_response")
        if not isinstance(response, dict):
            response = {}
        return {
            "identity": identity,
            "geometry_digest_valid": bool(geometry) and geometry_digest == digest(geometry),
            "frame_valid": frame_ok(frame),
            "frame_digest_valid": frame_ok(frame) and frame_digest == digest(frame),
            "geometry_digest": geometry_digest,
            "frame_digest": frame_digest,
            "refinement": refinement_signature(identity),
            "response": response,
        }

    candidate = lane_identity("candidate")
    reference = lane_identity("reference")

    identity_checks = {
        "canonical_geometry_digest_matches": candidate["geometry_digest_valid"]
        and reference["geometry_digest_valid"]
        and candidate["geometry_digest"] == reference["geometry_digest"],
        "right_handed_frame_digest_matches": candidate["frame_digest_valid"]
        and reference["frame_digest_valid"]
        and candidate["frame_digest"] == reference["frame_digest"],
        "unique_monotone_refinement_identity_matches": candidate["refinement"][0]
        and reference["refinement"][0]
        and candidate["refinement"][1] == reference["refinement"][1],
    }

    def derivative_residual(response: dict[str, Any]) -> dict[str, Any]:
        try:
            h_values = [float(value) for value in response.get("H_A_per_m")]
            b_values = [float(value) for value in response.get("B_T")]
            energy = [float(value) for value in response.get("energy_density_J_per_m3")]
            coenergy = [
                float(value) for value in response.get("coenergy_density_J_per_m3")
            ]
        except (TypeError, ValueError):
            return {"valid": False, "max_relative_residual": math.inf}
        valid_arrays = (
            len(h_values) >= sample_limit
            and len(h_values) == len(b_values) == len(energy) == len(coenergy)
            and all(
                math.isfinite(value)
                for values in (h_values, b_values, energy, coenergy)
                for value in values
            )
            and all(right > left for left, right in zip(h_values, h_values[1:]))
            and all(right > left for left, right in zip(b_values, b_values[1:]))
        )
        if not valid_arrays:
            return {"valid": False, "max_relative_residual": math.inf}
        residuals: list[float] = []
        for index in range(len(h_values) - 1):
            d_energy_d_b = (energy[index + 1] - energy[index]) / (
                b_values[index + 1] - b_values[index]
            )
            d_coenergy_d_h = (coenergy[index + 1] - coenergy[index]) / (
                h_values[index + 1] - h_values[index]
            )
            residuals.extend(
                [
                    abs(d_energy_d_b - 0.5 * (h_values[index] + h_values[index + 1]))
                    / max(abs(d_energy_d_b), 1.0),
                    abs(d_coenergy_d_h - 0.5 * (b_values[index] + b_values[index + 1]))
                    / max(abs(d_coenergy_d_h), 1.0),
                ]
            )
        maximum = max(residuals, default=math.inf)
        return {
            "valid": maximum <= residual_limit,
            "max_relative_residual": maximum,
            "sample_count": len(h_values),
        }

    candidate_derivative = derivative_residual(candidate["response"])
    reference_derivative = derivative_residual(reference["response"])
    checks = {
        **identity_checks,
        "candidate_energy_derivative_identity": candidate_derivative["valid"],
        "reference_energy_derivative_identity": reference_derivative["valid"],
    }
    return {
        "policy": "nonlinear_field_energy_identity_gate_v5",
        "status": "ok" if all(checks.values()) else "needs_attention",
        "accepted": all(checks.values()),
        "checks": checks,
        "issues": [name for name, accepted in checks.items() if not accepted],
        "identity": {
            "candidate_geometry_digest": candidate["geometry_digest"],
            "reference_geometry_digest": reference["geometry_digest"],
            "candidate_frame_digest": candidate["frame_digest"],
            "reference_frame_digest": reference["frame_digest"],
        },
        "derivative_identity": {
            "candidate": candidate_derivative,
            "reference": reference_derivative,
        },
        "tolerances": {
            "max_energy_derivative_relative_residual": residual_limit,
            "min_response_samples": sample_limit,
        },
        "notes": [
            "canonical geometry and frame digests are recomputed from supplied structures",
            "refinement levels require unique identifiers, unique mesh identities, and increasing element counts",
            "dW/dB=H and dWstar/dH=B are checked before any cross-lane numeric claim",
        ],
    }


def nonlinear_field_energy_artifact_contract_gate_v6(
    summary: dict[str, Any],
    *,
    max_residual_norm: float = 1.0e-8,
    min_response_samples: int = 4,
) -> dict[str, Any]:
    """Validate completed field-energy result artifacts before comparison.

    The v5 gate proves geometric and constitutive identities. This independent
    gate checks the result envelope itself: schema, completion/convergence,
    SI observable units, finite ordered samples, and a recomputable response
    digest. It deliberately does not read paths or claim solver parity.
    """

    if not isinstance(summary, dict):
        raise ValueError("summary must be a mapping")
    residual_limit = float(max_residual_norm)
    sample_limit = int(min_response_samples)
    if not math.isfinite(residual_limit) or residual_limit < 0.0:
        raise ValueError("max_residual_norm must be finite and nonnegative")
    if sample_limit < 4:
        raise ValueError("min_response_samples must be at least 4")

    def digest(value: object) -> str:
        encoded = json.dumps(
            value, ensure_ascii=True, sort_keys=True, separators=(",", ":")
        ).encode("utf-8")
        return hashlib.sha256(encoded).hexdigest()

    def valid_sha(value: object) -> bool:
        text = str(value or "").strip().lower()
        return len(text) == 64 and all(
            character in "0123456789abcdef" for character in text
        )

    required_columns = [
        "sample_id",
        "H_A_per_m",
        "B_T",
        "energy_density_J_per_m3",
        "coenergy_density_J_per_m3",
    ]
    required_units = {
        "sample_id": "1",
        "H_A_per_m": "A/m",
        "B_T": "T",
        "energy_density_J_per_m3": "J/m^3",
        "coenergy_density_J_per_m3": "J/m^3",
    }
    identity_keys = (
        "geometry_identity_sha256",
        "material_table_sha256",
        "excitation_identity_sha256",
        "coordinate_system",
        "unit_system",
    )

    def lane_contract(name: str) -> dict[str, Any]:
        raw = summary.get(name)
        if not isinstance(raw, dict):
            raise ValueError(f"{name} must be a mapping")
        identity = raw.get("identity")
        identity = identity if isinstance(identity, dict) else {}
        response = raw.get("response")
        response = response if isinstance(response, dict) else {}
        columns = response.get("observable_columns")
        units = response.get("observable_units")
        sample_ids = response.get("sample_id")
        h_values = response.get("H_A_per_m")
        b_values = response.get("B_T")
        energy = response.get("energy_density_J_per_m3")
        coenergy = response.get("coenergy_density_J_per_m3")
        arrays = (sample_ids, h_values, b_values, energy, coenergy)
        arrays_are_lists = all(isinstance(values, list) for values in arrays)
        lengths_match = arrays_are_lists and len({len(values) for values in arrays}) == 1
        sample_count = len(h_values) if isinstance(h_values, list) else 0
        finite = False
        ordered = False
        if lengths_match and sample_count >= sample_limit:
            try:
                finite = all(
                    math.isfinite(float(value))
                    for values in arrays[1:]
                    for value in values
                ) and all(isinstance(value, int) and not isinstance(value, bool) for value in sample_ids)
                ordered = (
                    len(set(sample_ids)) == len(sample_ids)
                    and all(right > left for left, right in zip(sample_ids, sample_ids[1:]))
                    and all(right > left for left, right in zip(h_values, h_values[1:]))
                    and all(right >= left for left, right in zip(b_values, b_values[1:]))
                )
            except (TypeError, ValueError):
                finite = False
                ordered = False
        response_core = {
            "observable_columns": columns,
            "observable_units": units,
            "dimension_order": response.get("dimension_order"),
            "sample_id": sample_ids,
            "H_A_per_m": h_values,
            "B_T": b_values,
            "energy_density_J_per_m3": energy,
            "coenergy_density_J_per_m3": coenergy,
        }
        checks = {
            "artifact_schema_supported": raw.get("schema") == "radia.nonlinear-field-energy-artifact.v1",
            "artifact_id_present": bool(str(raw.get("artifact_id") or "").strip()),
            "result_completed": raw.get("status") == "completed",
            "solver_converged": raw.get("solver_converged") is True,
            "residual_norm_finite_and_bounded": (
                isinstance(raw.get("residual_norm"), (int, float))
                and math.isfinite(float(raw["residual_norm"]))
                and float(raw["residual_norm"]) <= residual_limit
            ),
            "identity_contract_present": all(
                bool(str(identity.get(key) or "").strip()) for key in identity_keys
            ),
            "identity_digests_valid": all(
                valid_sha(identity.get(key))
                for key in identity_keys[:3]
            ),
            "observable_columns_exact": columns == required_columns,
            "observable_units_exact": units == required_units,
            "dimension_order_is_c": response.get("dimension_order") == "C",
            "response_arrays_are_lists": arrays_are_lists,
            "response_lengths_match": lengths_match,
            "response_sample_count_sufficient": sample_count >= sample_limit,
            "response_values_finite": finite,
            "response_samples_strict_and_monotone": ordered,
            "response_digest_matches": (
                valid_sha(identity.get("response_identity_sha256"))
                and identity.get("response_identity_sha256", "").casefold() == digest(response_core)
            ),
        }
        return {
            "checks": checks,
            "response_digest": str(identity.get("response_identity_sha256") or "").lower(),
            "sample_count": sample_count,
        }

    candidate = lane_contract("candidate")
    reference = lane_contract("reference")
    candidate_identity = summary["candidate"].get("identity", {})
    reference_identity = summary["reference"].get("identity", {})
    identity_matches = {
        key: candidate_identity.get(key) == reference_identity.get(key)
        for key in identity_keys
    }
    checks = {
        "candidate_contract_valid": all(candidate["checks"].values()),
        "reference_contract_valid": all(reference["checks"].values()),
        "artifact_ids_are_distinct": (
            summary["candidate"].get("artifact_id") != summary["reference"].get("artifact_id")
        ),
        "comparison_identity_matches": all(identity_matches.values()),
    }
    checks.update({f"comparison_{key}": value for key, value in identity_matches.items()})
    return {
        "policy": "nonlinear_field_energy_artifact_contract_gate_v6",
        "status": "ok" if all(checks.values()) else "needs_attention",
        "accepted": all(checks.values()),
        "checks": checks,
        "issues": [name for name, accepted in checks.items() if not accepted],
        "lane_details": {"candidate": candidate, "reference": reference},
        "tolerances": {
            "max_residual_norm": residual_limit,
            "min_response_samples": sample_limit,
        },
        "notes": [
            "result paths are never opened; only the supplied artifact contract is checked",
            "units and dimension order are explicit to prevent numerically plausible misinterpretation",
            "a valid result envelope is necessary but insufficient for cross-solver numerical parity",
        ],
    }


def nonlinear_field_energy_lineage_gate_v7(summary: dict[str, Any]) -> dict[str, Any]:
    """Validate result lineage and run identity before cross-solver comparison.

    This gate is independent of the v6 result-envelope checks. It verifies that
    each result is a solver output with a declared parent, that its run digest
    is recomputable from the physical identity, and that both lanes share the
    same geometry/material/excitation/mesh lineage.
    """

    if not isinstance(summary, dict):
        raise ValueError("summary must be a mapping")

    def digest(value: object) -> str:
        encoded = json.dumps(
            value, ensure_ascii=True, sort_keys=True, separators=(",", ":")
        ).encode("utf-8")
        return hashlib.sha256(encoded).hexdigest()

    def valid_sha(value: object) -> bool:
        text = str(value or "").strip().lower()
        return len(text) == 64 and all(
            character in "0123456789abcdef" for character in text
        )

    identity_keys = (
        "geometry_identity_sha256",
        "material_table_sha256",
        "excitation_identity_sha256",
        "mesh_identity_sha256",
    )
    run_keys = (*identity_keys, "solver_version", "run_id")

    def lane_contract(name: str) -> dict[str, Any]:
        lane = summary.get(name)
        if not isinstance(lane, dict):
            raise ValueError(f"{name} must be a mapping")
        identity = lane.get("identity")
        identity = identity if isinstance(identity, dict) else {}
        lineage = lane.get("lineage")
        lineage = lineage if isinstance(lineage, dict) else {}
        run_core = {key: identity.get(key) for key in run_keys}
        checks = {
            "artifact_id_present": bool(str(lane.get("artifact_id") or "").strip()),
            "lineage_is_solver_output": lineage.get("source_kind") == "solver_output",
            "root_case_id_present": bool(str(lineage.get("root_case_id") or "").strip()),
            "parent_artifact_id_present": bool(str(lineage.get("parent_artifact_id") or "").strip()),
            "parent_artifact_digest_valid": valid_sha(lineage.get("parent_artifact_sha256")),
            "identity_fields_present": all(
                bool(str(identity.get(key) or "").strip()) for key in run_keys
            ),
            "identity_digests_valid": all(valid_sha(identity.get(key)) for key in identity_keys),
            "run_identity_digest_matches": (
                valid_sha(identity.get("run_identity_sha256"))
                and str(identity.get("run_identity_sha256")).casefold() == digest(run_core)
            ),
        }
        return {
            "checks": checks,
            "root_case_id": str(lineage.get("root_case_id") or ""),
            "parent_artifact_id": str(lineage.get("parent_artifact_id") or ""),
            "run_id": str(identity.get("run_id") or ""),
        }

    candidate = lane_contract("candidate")
    reference = lane_contract("reference")
    candidate_lane = summary["candidate"]
    reference_lane = summary["reference"]
    candidate_identity = candidate_lane.get("identity", {})
    reference_identity = reference_lane.get("identity", {})
    candidate_lineage = candidate_lane.get("lineage", {})
    reference_lineage = reference_lane.get("lineage", {})
    identity_matches = {
        key: candidate_identity.get(key) == reference_identity.get(key)
        for key in identity_keys
    }
    checks = {
        "candidate_lineage_contract_valid": all(candidate["checks"].values()),
        "reference_lineage_contract_valid": all(reference["checks"].values()),
        "artifact_ids_are_distinct": candidate_lane.get("artifact_id") != reference_lane.get("artifact_id"),
        "run_ids_are_distinct": candidate["run_id"] != reference["run_id"],
        "root_case_identity_matches": (
            candidate["root_case_id"]
            and candidate["root_case_id"] == reference["root_case_id"]
        ),
        "parent_lineage_is_declared": (
            candidate_lineage.get("parent_artifact_id")
            and reference_lineage.get("parent_artifact_id")
        ),
        "cross_lane_physical_identity_matches": all(identity_matches.values()),
    }
    checks.update({f"cross_lane_{key}": value for key, value in identity_matches.items()})
    return {
        "policy": "nonlinear_field_energy_lineage_gate_v7",
        "status": "ok" if all(checks.values()) else "needs_attention",
        "accepted": all(checks.values()),
        "checks": checks,
        "issues": [name for name, accepted in checks.items() if not accepted],
        "lane_details": {"candidate": candidate, "reference": reference},
        "notes": [
            "run identity is recomputed from physical identity and solver run metadata",
            "parent lineage is metadata-only and does not open or trust arbitrary paths",
            "matching lineage is necessary but insufficient for numerical solver parity",
        ],
    }


def nonlinear_field_energy_physical_admissibility_gate_v8(
    summary: dict[str, Any],
    *,
    max_legendre_relative_residual: float = 1.0e-8,
    min_response_samples: int = 4,
) -> dict[str, Any]:
    """Reject physically inadmissible nonlinear field-energy responses."""

    if not isinstance(summary, dict):
        raise ValueError("summary must be a mapping")
    residual_limit = float(max_legendre_relative_residual)
    sample_limit = int(min_response_samples)
    if not math.isfinite(residual_limit) or residual_limit < 0.0:
        raise ValueError("max_legendre_relative_residual must be finite and nonnegative")
    if sample_limit < 4:
        raise ValueError("min_response_samples must be at least 4")

    def lane_contract(name: str) -> dict[str, Any]:
        lane = summary.get(name)
        if not isinstance(lane, dict):
            raise ValueError(f"{name} must be a mapping")
        response = lane.get("response")
        response = response if isinstance(response, dict) else {}
        names = (
            "H_A_per_m",
            "B_T",
            "energy_density_J_per_m3",
            "coenergy_density_J_per_m3",
            "differential_permeability_H_per_m",
        )
        arrays = {key: response.get(key) for key in names}
        arrays_are_lists = all(isinstance(values, list) for values in arrays.values())
        lengths_match = arrays_are_lists and len({len(values) for values in arrays.values()}) == 1
        finite = False
        ordered = False
        nonnegative = False
        positive_tangent = False
        legendre_residual = math.inf
        if lengths_match and len(arrays["H_A_per_m"]) >= sample_limit:
            try:
                finite = all(
                    math.isfinite(float(value))
                    for values in arrays.values()
                    for value in values
                )
                ordered = (
                    all(right > left for left, right in zip(arrays["H_A_per_m"], arrays["H_A_per_m"][1:]))
                    and all(right >= left for left, right in zip(arrays["B_T"], arrays["B_T"][1:]))
                )
                nonnegative = all(
                    float(value) >= 0.0
                    for key in ("energy_density_J_per_m3", "coenergy_density_J_per_m3")
                    for value in arrays[key]
                )
                positive_tangent = all(float(value) > 0.0 for value in arrays["differential_permeability_H_per_m"])
                residuals = [
                    abs(float(energy) + float(coenergy) - float(h) * float(b))
                    / max(abs(float(h) * float(b)), 1.0)
                    for h, b, energy, coenergy in zip(
                        arrays["H_A_per_m"],
                        arrays["B_T"],
                        arrays["energy_density_J_per_m3"],
                        arrays["coenergy_density_J_per_m3"],
                    )
                ]
                legendre_residual = max(residuals, default=math.inf)
            except (TypeError, ValueError, ZeroDivisionError):
                finite = ordered = nonnegative = positive_tangent = False
        checks = {
            "response_arrays_are_lists": arrays_are_lists,
            "response_lengths_match": lengths_match,
            "response_sample_count_sufficient": (
                isinstance(arrays["H_A_per_m"], list)
                and len(arrays["H_A_per_m"]) >= sample_limit
            ),
            "response_values_finite": finite,
            "field_samples_ordered": ordered,
            "energy_and_coenergy_nonnegative": nonnegative,
            "differential_permeability_positive": positive_tangent,
            "legendre_energy_identity_satisfied": legendre_residual <= residual_limit,
        }
        return {"checks": checks, "max_legendre_relative_residual": legendre_residual}

    candidate = lane_contract("candidate")
    reference = lane_contract("reference")
    candidate_case = summary["candidate"].get("identity", {}).get("comparison_case_id")
    reference_case = summary["reference"].get("identity", {}).get("comparison_case_id")
    checks = {
        "candidate_physical_contract_valid": all(candidate["checks"].values()),
        "reference_physical_contract_valid": all(reference["checks"].values()),
        "comparison_case_identity_matches": bool(candidate_case) and candidate_case == reference_case,
    }
    return {
        "policy": "nonlinear_field_energy_physical_admissibility_gate_v8",
        "status": "ok" if all(checks.values()) else "needs_attention",
        "accepted": all(checks.values()),
        "checks": checks,
        "issues": [name for name, accepted in checks.items() if not accepted],
        "lane_details": {"candidate": candidate, "reference": reference},
        "tolerances": {
            "max_legendre_relative_residual": residual_limit,
            "min_response_samples": sample_limit,
        },
        "notes": [
            "energy and coenergy must be nonnegative with positive differential permeability",
            "W + W* = H dot B is recomputed at every response sample",
            "physical admissibility is necessary but insufficient for cross-solver numerical parity",
        ],
    }


def nonlinear_field_energy_observable_comparison_gate_v9(
    summary: dict[str, Any],
    *,
    max_relative_difference: float = 0.05,
    absolute_floor: float = 1.0e-12,
) -> dict[str, Any]:
    """Compare observables only after identity and sample contracts match."""

    if not isinstance(summary, dict):
        raise ValueError("summary must be a mapping")
    difference_limit = float(max_relative_difference)
    floor = float(absolute_floor)
    if not math.isfinite(difference_limit) or difference_limit < 0.0:
        raise ValueError("max_relative_difference must be finite and nonnegative")
    if not math.isfinite(floor) or floor <= 0.0:
        raise ValueError("absolute_floor must be finite and positive")
    identity_keys = (
        "geometry_identity_sha256",
        "material_table_sha256",
        "excitation_identity_sha256",
        "mesh_identity_sha256",
        "coordinate_system",
        "unit_system",
    )
    observable_keys = ("average_B_T", "rms_B_T", "energy_J", "coenergy_J")
    required_units = {"average_B_T": "T", "rms_B_T": "T", "energy_J": "J", "coenergy_J": "J"}

    def valid_sha(value: object) -> bool:
        text = str(value or "").strip().lower()
        return len(text) == 64 and all(character in "0123456789abcdef" for character in text)

    def lane_contract(name: str) -> dict[str, Any]:
        lane = summary.get(name)
        if not isinstance(lane, dict):
            raise ValueError(f"{name} must be a mapping")
        identity = lane.get("identity")
        identity = identity if isinstance(identity, dict) else {}
        observables = lane.get("observables")
        observables = observables if isinstance(observables, dict) else {}
        sample_ids = observables.get("sample_id")
        arrays = {key: observables.get(key) for key in observable_keys}
        arrays_are_lists = isinstance(sample_ids, list) and all(isinstance(value, list) for value in arrays.values())
        lengths_match = arrays_are_lists and len({len(sample_ids), *(len(value) for value in arrays.values())}) == 1
        finite = False
        sample_contract = False
        max_difference = math.inf
        if lengths_match:
            try:
                finite = all(math.isfinite(float(item)) for values in arrays.values() for item in values)
                sample_contract = (
                    bool(sample_ids)
                    and all(isinstance(item, (int, str)) and not isinstance(item, bool) for item in sample_ids)
                    and len(set(sample_ids)) == len(sample_ids)
                )
            except (TypeError, ValueError):
                finite = sample_contract = False
        units = observables.get("observable_units")
        checks = {
            "identity_fields_present": all(bool(str(identity.get(key) or "").strip()) for key in identity_keys),
            "identity_digests_valid": all(valid_sha(identity.get(key)) for key in identity_keys[:4]),
            "observable_columns_exact": observables.get("observable_columns") == list(observable_keys),
            "observable_units_exact": units == required_units,
            "sample_arrays_are_lists": arrays_are_lists,
            "sample_array_lengths_match": lengths_match,
            "observable_values_finite": finite,
            "sample_identity_valid": sample_contract,
        }
        return {"identity": identity, "observables": observables, "checks": checks, "sample_ids": sample_ids, "arrays": arrays, "max_relative_difference": max_difference}

    candidate = lane_contract("candidate")
    reference = lane_contract("reference")
    identity_matches = {key: candidate["identity"].get(key) == reference["identity"].get(key) for key in identity_keys}
    sample_ids_match = candidate["sample_ids"] == reference["sample_ids"]
    if (
        candidate["checks"]["observable_values_finite"]
        and reference["checks"]["observable_values_finite"]
        and candidate["sample_ids"]
        and reference["sample_ids"]
        and sample_ids_match
    ):
        differences = [
            abs(float(left) - float(right)) / max(abs(float(left)), abs(float(right)), floor)
            for key in observable_keys
            for left, right in zip(candidate["arrays"][key], reference["arrays"][key])
        ]
        max_difference = max(differences, default=math.inf)
    else:
        max_difference = math.inf
    checks = {
        "candidate_contract_valid": all(candidate["checks"].values()),
        "reference_contract_valid": all(reference["checks"].values()),
        "comparison_identity_matches": all(identity_matches.values()),
        "observable_sample_identity_matches": sample_ids_match,
        "observable_differences_within_limit": max_difference <= difference_limit,
    }
    return {
        "policy": "nonlinear_field_energy_observable_comparison_gate_v9",
        "status": "ok" if all(checks.values()) else "needs_attention",
        "accepted": all(checks.values()),
        "checks": checks,
        "issues": [name for name, accepted in checks.items() if not accepted],
        "max_relative_difference": max_difference,
        "tolerances": {"max_relative_difference": difference_limit, "absolute_floor": floor},
        "notes": [
            "identity and sample contracts are checked before numeric differences",
            "relative differences use a positive absolute floor for near-zero observables",
            "numeric agreement is not a claim of solver correctness without independent physics gates",
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


def nonlinear_constitutive_point_sample_gate(
    summary: dict[str, Any],
    *,
    max_response_relative_difference: float = 0.05,
    max_direction_sine: float = 0.02,
    max_vector_magnitude_relative_residual: float = 1.0e-9,
    min_response_samples: int = 5,
) -> dict[str, Any]:
    """Gate pointwise constitutive evidence without accepting recovered fields.

    A postprocessor's nodally averaged or equivalent-source-integrated B/H
    values are useful spatial-field evidence, but they are not an unsmoothed
    element-local evaluation of the material law.  This gate keeps that
    distinction explicit and only admits an identity-bound ``ELEMENT_LOCAL``
    reference as constitutive evidence.
    """

    if not isinstance(summary, dict):
        raise ValueError("summary must be a mapping")
    response_limit = float(max_response_relative_difference)
    direction_limit = float(max_direction_sine)
    magnitude_limit = float(max_vector_magnitude_relative_residual)
    sample_limit = int(min_response_samples)
    if any(
        not math.isfinite(value) or value < 0.0
        for value in (response_limit, direction_limit, magnitude_limit)
    ):
        raise ValueError("relative tolerances must be finite and nonnegative")
    if direction_limit > 1.0:
        raise ValueError("max_direction_sine must not exceed 1")
    if sample_limit < 5:
        raise ValueError("min_response_samples must be at least 5")

    source = summary.get("source")
    candidate = summary.get("candidate")
    source = source if isinstance(source, dict) else {}
    candidate = candidate if isinstance(candidate, dict) else {}
    source_identity = source.get("identity")
    source_identity = source_identity if isinstance(source_identity, dict) else {}
    candidate_identity = candidate.get("identity")
    candidate_identity = (
        candidate_identity if isinstance(candidate_identity, dict) else {}
    )
    values = source.get("values")
    values = values if isinstance(values, list) else []

    def canonical_digest(value: Any) -> str:
        payload = json.dumps(
            value, ensure_ascii=True, separators=(",", ":"), sort_keys=True
        ).encode("ascii")
        return hashlib.sha256(payload).hexdigest()

    source_rows: list[tuple[float, float, float, float]] = []
    source_schema_valid = bool(values)
    max_b_magnitude_residual = math.inf
    max_h_magnitude_residual = math.inf
    max_cross_relative = math.inf
    min_cosine = -1.0
    b_residuals: list[float] = []
    h_residuals: list[float] = []
    cross_values: list[float] = []
    cosine_values: list[float] = []
    for value in values:
        if not isinstance(value, dict):
            source_schema_valid = False
            break
        try:
            b_vector = [float(item) for item in value["field_T"]]
            h_vector = [
                float(item)
                for item in value["magnetic_field_strength_A_per_m"]
            ]
            b_magnitude = float(value["magnitude_T"])
            h_magnitude = float(
                value["magnetic_field_strength_magnitude_A_per_m"]
            )
        except (KeyError, TypeError, ValueError):
            source_schema_valid = False
            break
        if len(b_vector) != 3 or len(h_vector) != 3 or not all(
            math.isfinite(item)
            for item in (*b_vector, *h_vector, b_magnitude, h_magnitude)
        ):
            source_schema_valid = False
            break
        b_norm = math.sqrt(sum(item * item for item in b_vector))
        h_norm = math.sqrt(sum(item * item for item in h_vector))
        denominator = b_norm * h_norm
        if denominator <= 1.0e-300:
            source_schema_valid = False
            break
        cross = (
            (b_vector[1] * h_vector[2] - b_vector[2] * h_vector[1]) ** 2
            + (b_vector[2] * h_vector[0] - b_vector[0] * h_vector[2]) ** 2
            + (b_vector[0] * h_vector[1] - b_vector[1] * h_vector[0]) ** 2
        ) ** 0.5 / denominator
        cosine = sum(b * h for b, h in zip(b_vector, h_vector)) / denominator
        b_residuals.append(abs(b_norm - b_magnitude) / max(abs(b_magnitude), 1.0e-300))
        h_residuals.append(abs(h_norm - h_magnitude) / max(abs(h_magnitude), 1.0e-300))
        cross_values.append(cross)
        cosine_values.append(cosine)
        source_rows.append((h_magnitude, b_magnitude, cross, cosine))

    if source_schema_valid:
        source_rows.sort(key=lambda row: row[0])
        max_b_magnitude_residual = max(b_residuals)
        max_h_magnitude_residual = max(h_residuals)
        max_cross_relative = max(cross_values)
        min_cosine = min(cosine_values)

    try:
        candidate_h = [float(item) for item in candidate.get("H_A_per_m", [])]
        candidate_b = [float(item) for item in candidate.get("B_T", [])]
    except (TypeError, ValueError):
        candidate_h = []
        candidate_b = []
    candidate_valid = (
        len(candidate_h) == len(candidate_b)
        and len(candidate_h) >= sample_limit
        and all(math.isfinite(item) for item in (*candidate_h, *candidate_b))
        and all(right > left for left, right in zip(candidate_h, candidate_h[1:]))
    )
    source_h = [row[0] for row in source_rows]
    source_b = [row[1] for row in source_rows]
    candidate_at_source: list[float] = []
    grids_match = candidate_valid
    if grids_match:
        for source_value in source_h:
            matching = [
                b_value
                for h_value, b_value in zip(candidate_h, candidate_b)
                if math.isclose(
                    source_value, h_value, rel_tol=1.0e-12, abs_tol=1.0e-12
                )
            ]
            if len(matching) != 1:
                grids_match = False
                candidate_at_source = []
                break
            candidate_at_source.append(matching[0])
    pointwise_relative_difference: list[float] = []
    rms_relative_difference: float | None = None
    if source_schema_valid and grids_match:
        pointwise_relative_difference = [
            abs(actual - expected) / max(abs(expected), 1.0e-300)
            for actual, expected in zip(source_b, candidate_at_source)
        ]
        rms_relative_difference = math.sqrt(
            sum(
                (actual - expected) ** 2
                for actual, expected in zip(source_b, candidate_at_source)
            )
            / max(
                sum(expected * expected for expected in candidate_at_source),
                1.0e-300,
            )
        )

    point_set = source_identity.get("point_set_m")
    point_set = point_set if isinstance(point_set, list) else []
    material_curve = source_identity.get("material_curve")
    material_curve = material_curve if isinstance(material_curve, dict) else {}
    claimed_identity_digest = str(
        source_identity.get("constitutive_identity_sha256") or ""
    )
    identity_core = {
        key: value
        for key, value in source_identity.items()
        if key != "constitutive_identity_sha256"
    }
    checks = {
        "source_completed": source.get("status") == "completed",
        "sample_count_sufficient": len(source_rows) >= sample_limit,
        "joint_point_fields_finite": source_schema_valid,
        "SI_units_explicit": source_identity.get("response_units")
        == {"H": "A/m", "B": "T"},
        "point_set_identity_valid": bool(point_set)
        and source_identity.get("point_set_sha256") == canonical_digest(point_set),
        "source_identity_digest_valid": bool(identity_core)
        and claimed_identity_digest == canonical_digest(identity_core),
        "material_table_identity_matches": bool(
            material_curve.get("canonical_table_sha256")
        )
        and material_curve.get("canonical_table_sha256")
        == candidate_identity.get("bh_table_sha256"),
        "source_is_unsmoothed_element_local": (
            source.get("response_evidence_status")
            == "source_native_element_local_B_H_samples"
            and source.get("constitutive_oracle") is True
            and source.get("constitutive_comparison_ready") is True
            and source_identity.get("constitutive_oracle") is True
            and source_identity.get("field_recovery") == "ELEMENT_LOCAL"
            and source_identity.get("field_recovery_semantics")
            == "unsmoothed_element_local_evaluation"
        ),
        "candidate_grid_matches_source_samples": grids_match,
        "reported_vector_magnitudes_match": (
            max_b_magnitude_residual <= magnitude_limit
            and max_h_magnitude_residual <= magnitude_limit
        ),
        "isotropic_B_H_directions_collinear": (
            max_cross_relative <= direction_limit and min_cosine >= 0.0
        ),
        "B_response_matches": bool(pointwise_relative_difference)
        and max(pointwise_relative_difference) <= response_limit
        and rms_relative_difference is not None
        and rms_relative_difference <= response_limit,
    }
    return {
        "policy": "nonlinear_constitutive_point_sample_gate_v1",
        "status": "ok" if all(checks.values()) else "needs_attention",
        "checks": checks,
        "issues": [name for name, accepted in checks.items() if not accepted],
        "comparison_performed": source_schema_valid and grids_match,
        "constitutive_parity_established": all(checks.values()),
        "metrics": {
            "sample_count": len(source_rows),
            "H_range_A_per_m": (
                [min(source_h), max(source_h)] if source_h else None
            ),
            "max_B_relative_difference": (
                max(pointwise_relative_difference)
                if pointwise_relative_difference
                else None
            ),
            "rms_B_relative_difference": rms_relative_difference,
            "max_B_H_direction_sine": max_cross_relative,
            "min_B_H_direction_cosine": min_cosine,
            "max_B_magnitude_relative_residual": max_b_magnitude_residual,
            "max_H_magnitude_relative_residual": max_h_magnitude_residual,
        },
        "source_recovery": {
            "mode": source_identity.get("field_recovery"),
            "semantics": source_identity.get("field_recovery_semantics"),
            "constitutive_oracle": source_identity.get("constitutive_oracle"),
        },
        "tolerances": {
            "response_relative_difference": response_limit,
            "direction_sine": direction_limit,
            "vector_magnitude_relative_residual": magnitude_limit,
            "minimum_response_samples": sample_limit,
        },
        "notes": [
            "NODAL and INTEGRATION point fields remain valid spatial observables but are rejected as element-local material-law evidence",
            "an accepted point contract proves only the sampled constitutive response, not whole-solver equivalence",
        ],
    }


def controlled_uniform_field_constitutive_sweep_gate(
    summary: dict[str, Any],
    *,
    max_response_relative_difference: float = 0.01,
    min_response_samples: int = 5,
) -> dict[str, Any]:
    """Compare a bounded homogeneous-interior sweep with one B-H response.

    The source input must already have passed an independent multicase gate
    that binds every case to one material table, excitation, result database,
    and material region. This admits recovered NODAL fields only as a
    controlled experiment; it never promotes them to a general element-local
    constitutive oracle.
    """

    if not isinstance(summary, dict):
        raise ValueError("summary must be a mapping")
    response_limit = float(max_response_relative_difference)
    sample_limit = int(min_response_samples)
    if not math.isfinite(response_limit) or response_limit < 0.0:
        raise ValueError(
            "max_response_relative_difference must be finite and nonnegative"
        )
    if sample_limit < 5:
        raise ValueError("min_response_samples must be at least 5")

    source = summary.get("source_control")
    source = source if isinstance(source, dict) else {}
    candidate = summary.get("candidate")
    candidate = candidate if isinstance(candidate, dict) else {}
    source_identity = source.get("identity")
    source_identity = source_identity if isinstance(source_identity, dict) else {}
    candidate_identity = candidate.get("identity")
    candidate_identity = (
        candidate_identity if isinstance(candidate_identity, dict) else {}
    )
    cases = source.get("case_summaries")
    cases = cases if isinstance(cases, list) else []

    def is_sha256(value: Any) -> bool:
        text = str(value or "").lower()
        return len(text) == 64 and all(character in "0123456789abcdef" for character in text)

    parsed_cases: list[tuple[float, float]] = []
    case_schema_valid = bool(cases)
    for case in cases:
        if not isinstance(case, dict) or case.get("accepted") is not True:
            case_schema_valid = False
            break
        try:
            h_value = float(case["mean_H_A_per_m"])
            b_value = float(case["mean_B_T"])
        except (KeyError, TypeError, ValueError):
            case_schema_valid = False
            break
        if not math.isfinite(h_value) or not math.isfinite(b_value) or h_value <= 0.0:
            case_schema_valid = False
            break
        parsed_cases.append((h_value, b_value))

    parsed_cases.sort()
    source_h: list[float] = []
    source_b: list[float] = []
    duplicate_cases_consistent = True
    for h_value, b_value in parsed_cases:
        if source_h and math.isclose(
            h_value, source_h[-1], rel_tol=1.0e-12, abs_tol=1.0e-12
        ):
            duplicate_cases_consistent = duplicate_cases_consistent and math.isclose(
                b_value, source_b[-1], rel_tol=1.0e-9, abs_tol=1.0e-12
            )
            continue
        source_h.append(h_value)
        source_b.append(b_value)

    try:
        candidate_h = [float(value) for value in candidate.get("H_A_per_m", [])]
        candidate_b = [float(value) for value in candidate.get("B_T", [])]
    except (TypeError, ValueError):
        candidate_h = []
        candidate_b = []
    candidate_valid = (
        len(candidate_h) == len(candidate_b)
        and len(candidate_h) >= sample_limit
        and all(math.isfinite(value) for value in (*candidate_h, *candidate_b))
        and all(right > left for left, right in zip(candidate_h, candidate_h[1:]))
        and all(right >= left for left, right in zip(candidate_b, candidate_b[1:]))
    )
    grids_match = (
        candidate_valid
        and len(source_h) == len(candidate_h)
        and all(
            math.isclose(actual, expected, rel_tol=1.0e-12, abs_tol=1.0e-12)
            for actual, expected in zip(source_h, candidate_h)
        )
    )
    pointwise_relative_difference: list[float] = []
    rms_relative_difference: float | None = None
    if case_schema_valid and duplicate_cases_consistent and grids_match:
        pointwise_relative_difference = [
            abs(actual - expected) / max(abs(expected), 1.0e-300)
            for actual, expected in zip(source_b, candidate_b)
        ]
        rms_relative_difference = math.sqrt(
            sum((actual - expected) ** 2 for actual, expected in zip(source_b, candidate_b))
            / max(sum(expected * expected for expected in candidate_b), 1.0e-300)
        )

    required_source_digests = (
        source_identity.get("canonical_table_sha256"),
        source_identity.get("excitation_identity_sha256"),
        source_identity.get("result_database_sha256"),
        source_identity.get("case_set_sha256"),
    )
    source_checks = source.get("checks")
    source_checks = source_checks if isinstance(source_checks, dict) else {}
    source_metrics = source.get("metrics")
    source_metrics = source_metrics if isinstance(source_metrics, dict) else {}
    try:
        h_span_ratio = float(source_metrics.get("H_span_ratio", 0.0))
    except (TypeError, ValueError):
        h_span_ratio = 0.0
    checks = {
        "source_control_accepted": (
            source.get("status") == "accepted_as_constitutive_control"
            and source.get("constitutive_control_ready") is True
            and source_checks
            and all(value is True for value in source_checks.values())
        ),
        "source_claim_boundary_preserved": source.get("constitutive_oracle") is False,
        "source_identity_complete": all(is_sha256(value) for value in required_source_digests)
        and isinstance(source_identity.get("region_labels"), list)
        and bool(source_identity.get("region_labels")),
        "source_case_schema_valid": case_schema_valid,
        "duplicate_cases_consistent": duplicate_cases_consistent,
        "source_nonlinear_range_sufficient": len(source_h) >= sample_limit
        and h_span_ratio >= 100.0,
        "candidate_response_valid": candidate_valid,
        "material_table_identity_matches": (
            is_sha256(candidate_identity.get("bh_table_sha256"))
            and source_identity.get("canonical_table_sha256")
            == candidate_identity.get("bh_table_sha256")
        ),
        "candidate_grid_matches_source_control": grids_match,
        "B_response_matches": bool(pointwise_relative_difference)
        and max(pointwise_relative_difference) <= response_limit,
    }
    return {
        "policy": "controlled_uniform_field_constitutive_sweep_gate_v1",
        "status": "ok" if all(checks.values()) else "needs_attention",
        "checks": checks,
        "issues": [name for name, accepted in checks.items() if not accepted],
        "comparison_performed": bool(pointwise_relative_difference),
        "constitutive_control_parity_established": all(checks.values()),
        "metrics": {
            "source_case_count": len(parsed_cases),
            "distinct_H_count": len(source_h),
            "H_span_ratio": h_span_ratio,
            "max_B_relative_difference": (
                max(pointwise_relative_difference)
                if pointwise_relative_difference
                else None
            ),
            "rms_B_relative_difference": rms_relative_difference,
        },
        "identity": {
            "canonical_table_sha256": source_identity.get(
                "canonical_table_sha256"
            ),
            "excitation_identity_sha256": source_identity.get(
                "excitation_identity_sha256"
            ),
            "result_database_sha256": source_identity.get(
                "result_database_sha256"
            ),
            "case_set_sha256": source_identity.get("case_set_sha256"),
            "candidate_interpolation": candidate_identity.get(
                "constitutive_interpolation"
            ),
            "candidate_extrapolation": candidate_identity.get(
                "constitutive_extrapolation"
            ),
        },
        "tolerances": {
            "response_relative_difference": response_limit,
            "minimum_response_samples": sample_limit,
        },
        "notes": [
            "acceptance applies only to the identity-bound homogeneous-interior multicase control",
            "recovered NODAL fields remain unsuitable as a general element-local material-law oracle",
            "a failed response check identifies interpolation or solve-response mismatch without assigning either implementation as ground truth",
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
