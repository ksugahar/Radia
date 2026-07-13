"""Validation gate for families of linear two-winding inductance matrices."""
from __future__ import annotations

import math


def _relative_error(actual: float, expected: float) -> float:
    return abs(actual - expected) / max(abs(actual), abs(expected), 1.0e-300)


def inductance_matrix_family_gate(
    cases,
    *,
    expected_strongest_coupling_case: str | None = None,
    max_reciprocity_relative_error: float = 0.02,
    psd_relative_tolerance: float = 1.0e-12,
    max_identity_relative_error: float = 1.0e-6,
    max_replay_relative_error: float = 1.0e-9,
    max_turn_scaling_relative_error: float = 0.02,
):
    if not isinstance(cases, list) or len(cases) < 2:
        raise ValueError("cases must contain at least two matrix rows")
    if any(
        tolerance < 0.0
        for tolerance in (
            max_reciprocity_relative_error,
            psd_relative_tolerance,
            max_identity_relative_error,
            max_replay_relative_error,
            max_turn_scaling_relative_error,
        )
    ):
        raise ValueError("tolerances must be nonnegative")

    rows = []
    ids = []
    enrichment_flags = []
    for index, case in enumerate(cases):
        if not isinstance(case, dict):
            raise ValueError(f"case {index} must be an object")
        case_id = str(case.get("case_id") or "").strip()
        matrix = case.get("matrix_H")
        if not case_id or not isinstance(matrix, list) or len(matrix) != 2:
            raise ValueError(f"case {index} needs case_id and a 2x2 matrix_H")
        if any(not isinstance(row, list) or len(row) != 2 for row in matrix):
            raise ValueError(f"case {case_id} matrix_H must be 2x2")
        l11, m12 = (float(value) for value in matrix[0])
        m21, l22 = (float(value) for value in matrix[1])
        finite = all(math.isfinite(value) for value in (l11, m12, m21, l22))
        mutual = 0.5 * (m12 + m21)
        reciprocity = abs(m12 - m21) / max(abs(m12), abs(m21), 1.0e-300)
        diagonal_product = l11 * l22
        determinant = diagonal_product - mutual * mutual
        coupling = (
            abs(mutual) / math.sqrt(diagonal_product)
            if finite and l11 > 0.0 and l22 > 0.0
            else math.inf
        )
        enriched_fields = (
            "turns",
            "current_A",
            "flux_linkage_Vs",
            "energy_J",
            "coenergy_J",
            "replay_count",
            "replay_max_relative_error",
        )
        present = [field in case for field in enriched_fields]
        if any(present) and not all(present):
            raise ValueError(
                f"case {case_id} must provide all enriched identity/replay fields"
            )
        enriched = all(present)
        enrichment_flags.append(enriched)
        checks = {
            "all_finite": finite,
            "positive_self_inductances": l11 > 0.0 and l22 > 0.0,
            "mutual_terms_have_consistent_sign": m12 == 0.0 or m21 == 0.0 or m12 * m21 > 0.0,
            "reciprocity_within_tolerance": reciprocity <= max_reciprocity_relative_error,
            "symmetrized_matrix_positive_semidefinite": determinant
            >= -psd_relative_tolerance * max(abs(diagonal_product), 1.0e-300),
            "coupling_coefficient_bounded": coupling <= 1.0 + psd_relative_tolerance,
        }
        identity_metrics = None
        turns = None
        if enriched:
            turns = case["turns"]
            currents = case["current_A"]
            flux = case["flux_linkage_Vs"]
            if (
                not isinstance(turns, list)
                or len(turns) != 2
                or not isinstance(currents, list)
                or len(currents) != 2
                or not isinstance(flux, list)
                or len(flux) != 2
            ):
                raise ValueError(
                    f"case {case_id} turns/current_A/flux_linkage_Vs must be length two"
                )
            turns = [float(value) for value in turns]
            i1, i2 = (float(value) for value in currents)
            psi1, psi2 = (float(value) for value in flux)
            energy = float(case["energy_J"])
            coenergy = float(case["coenergy_J"])
            replay_count = int(case["replay_count"])
            replay_error = float(case["replay_max_relative_error"])
            predicted_flux = [l11 * i1 + m12 * i2, m21 * i1 + l22 * i2]
            flux_error = max(
                _relative_error(psi1, predicted_flux[0]),
                _relative_error(psi2, predicted_flux[1]),
            )
            predicted_energy = 0.5 * (
                l11 * i1 * i1 + (m12 + m21) * i1 * i2 + l22 * i2 * i2
            )
            energy_error = _relative_error(energy, predicted_energy)
            energy_coenergy_error = _relative_error(energy, coenergy)
            finite_identity = all(
                math.isfinite(value)
                for value in (
                    *turns,
                    i1,
                    i2,
                    psi1,
                    psi2,
                    energy,
                    coenergy,
                    replay_error,
                )
            )
            checks.update(
                {
                    "identity_values_finite": finite_identity,
                    "turn_counts_positive": all(value > 0.0 for value in turns),
                    "flux_linkage_matches_matrix_current": flux_error
                    <= max_identity_relative_error,
                    "energy_matches_quadratic_form": energy_error
                    <= max_identity_relative_error,
                    "linear_energy_matches_coenergy": energy_coenergy_error
                    <= max_identity_relative_error,
                    "independent_replay_is_stable": replay_count >= 2
                    and replay_error <= max_replay_relative_error,
                }
            )
            identity_metrics = {
                "flux_linkage_relative_error": flux_error,
                "energy_relative_error": energy_error,
                "energy_coenergy_relative_error": energy_coenergy_error,
                "replay_count": replay_count,
                "replay_max_relative_error": replay_error,
            }
        ids.append(case_id)
        rows.append(
            {
                "case_id": case_id,
                "topology_class": case.get("topology_class"),
                "matrix_H": [[l11, m12], [m21, l22]],
                "mutual_mean_H": mutual,
                "reciprocity_relative_error": reciprocity,
                "determinant_H2": determinant,
                "coupling_abs": coupling,
                "turns": turns,
                "identity_metrics": identity_metrics,
                "checks": checks,
                "status": "ok" if all(checks.values()) else "needs_attention",
            }
        )

    strongest = max(rows, key=lambda row: row["coupling_abs"])
    unique_ids = len(set(ids)) == len(ids)
    strongest_ok = (
        expected_strongest_coupling_case is None
        or strongest["case_id"] == expected_strongest_coupling_case
    )
    if any(enrichment_flags) and not all(enrichment_flags):
        raise ValueError("either every case or no case must use enriched identity fields")
    enriched_family = all(enrichment_flags)
    turn_scaling_checks = {}
    if enriched_family:
        ordered = sorted(rows, key=lambda row: row["turns"][1])
        baseline = ordered[0]
        turn_scaling_checks["primary_turn_count_is_fixed"] = all(
            row["turns"][0] == baseline["turns"][0] for row in ordered
        )
        turn_scaling_checks["secondary_turn_counts_are_distinct"] = len(
            {row["turns"][1] for row in ordered}
        ) == len(ordered)
        scaling_errors = []
        for row in ordered[1:]:
            denominators = (
                baseline["turns"][1],
                baseline["matrix_H"][1][1],
                baseline["mutual_mean_H"],
                baseline["coupling_abs"],
            )
            if not all(math.isfinite(value) and abs(value) > 0.0 for value in denominators):
                scaling_errors.append(math.inf)
                continue
            ratio = row["turns"][1] / baseline["turns"][1]
            scaling_errors.extend(
                [
                    _relative_error(row["matrix_H"][0][0], baseline["matrix_H"][0][0]),
                    _relative_error(
                        row["matrix_H"][1][1] / baseline["matrix_H"][1][1],
                        ratio * ratio,
                    ),
                    _relative_error(
                        abs(row["mutual_mean_H"] / baseline["mutual_mean_H"]), ratio
                    ),
                    _relative_error(row["coupling_abs"], baseline["coupling_abs"]),
                ]
            )
        turn_scaling_checks["turn_scaling_within_tolerance"] = bool(scaling_errors) and max(
            scaling_errors
        ) <= max_turn_scaling_relative_error
    else:
        scaling_errors = []
    family_checks = {
        "case_ids_unique": unique_ids,
        "all_case_matrices_valid": all(row["status"] == "ok" for row in rows),
        "expected_case_has_strongest_coupling": strongest_ok,
        **turn_scaling_checks,
    }
    return {
        "policy": (
            "inductance_matrix_family_gate_v2"
            if enriched_family
            else "inductance_matrix_family_gate_v1"
        ),
        "status": "ok" if all(family_checks.values()) else "needs_attention",
        "case_count": len(rows),
        "strongest_coupling_case": strongest["case_id"],
        "strongest_coupling_abs": strongest["coupling_abs"],
        "maximum_reciprocity_relative_error": max(
            row["reciprocity_relative_error"] for row in rows
        ),
        "maximum_turn_scaling_relative_error": max(scaling_errors, default=None),
        "checks": family_checks,
        "cases": rows,
        "lesson": (
            "Build a two-winding inductance matrix by exciting each winding separately. "
            "Gate reciprocal mutual terms, positive self terms, positive semidefiniteness, "
            "and |k|<=1 before comparing magnetic-circuit topologies. For turn-ratio "
            "families, also close psi=L I, W=I^T L I/2, W=W', independent replay, "
            "and the N, N^2 inductance scaling laws."
        ),
    }
