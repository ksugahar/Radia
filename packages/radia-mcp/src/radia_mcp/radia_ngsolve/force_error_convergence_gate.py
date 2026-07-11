"""Solver-neutral convergence gate for force-error profiles from multiple formulations."""

from __future__ import annotations

import math
from typing import Any


def dual_formulation_force_error_convergence_gate(
    formulation_rows: list[dict[str, Any]],
    *,
    reference_force: float,
    max_final_relative_error: float = 0.02,
    min_initial_to_final_improvement: float = 1.1,
    max_final_to_best_error_ratio: float = 1.5,
    max_tail_relative_span: float = 0.005,
) -> dict[str, Any]:
    """Gate final force accuracy and convergence envelopes without requiring monotonicity."""

    reference = float(reference_force)
    tolerances = [
        float(max_final_relative_error),
        float(min_initial_to_final_improvement),
        float(max_final_to_best_error_ratio),
        float(max_tail_relative_span),
    ]
    if not math.isfinite(reference) or reference <= 0.0:
        raise ValueError("reference_force must be finite and positive")
    if len(formulation_rows) < 2:
        raise ValueError("at least two formulation rows are required")
    if any(not math.isfinite(value) or value < 0.0 for value in tolerances):
        raise ValueError("tolerances and ratios must be finite and nonnegative")

    normalized = []
    for row_index, row in enumerate(formulation_rows):
        try:
            identifier = str(row["id"]).strip()
            levels = [float(value) for value in row["refinement_levels"]]
            observable_ids = [str(value).strip() for value in row["observable_ids"]]
            profiles = [
                [float(value) for value in profile] for profile in row["error_profiles"]
            ]
        except (KeyError, TypeError, ValueError) as exc:
            raise ValueError(f"formulation row {row_index} is malformed") from exc
        if not identifier:
            raise ValueError(f"formulation row {row_index} has an empty id")
        if len(levels) < 3:
            raise ValueError(f"formulation {identifier} needs at least three levels")
        if len(observable_ids) != len(profiles) or not profiles:
            raise ValueError(f"formulation {identifier} observable/profile counts differ")
        if len(set(observable_ids)) != len(observable_ids) or any(not value for value in observable_ids):
            raise ValueError(f"formulation {identifier} observable ids must be nonempty and unique")
        if any(len(profile) != len(levels) for profile in profiles):
            raise ValueError(f"formulation {identifier} profile lengths must match levels")
        if not all(math.isfinite(value) for value in levels + [x for p in profiles for x in p]):
            raise ValueError(f"formulation {identifier} contains a non-finite value")
        normalized.append({
            "id": identifier,
            "refinement_levels": levels,
            "observable_ids": observable_ids,
            "error_profiles": profiles,
        })

    summaries = []
    for row in normalized:
        profile_summaries = []
        for observable_id, profile in zip(row["observable_ids"], row["error_profiles"]):
            best = min(profile)
            final = profile[-1]
            initial = profile[0]
            improvement = initial / max(final, math.ulp(1.0))
            final_to_best = final / max(best, math.ulp(1.0))
            reversals = sum(
                right > left for left, right in zip(profile, profile[1:])
            )
            tail_span_relative = (
                abs(profile[-1] - profile[-2]) / reference
                if len(profile) >= 4
                else None
            )
            profile_summaries.append({
                "observable_id": observable_id,
                "initial_error": initial,
                "final_error": final,
                "best_error": best,
                "final_relative_error": final / reference,
                "initial_to_final_improvement": improvement,
                "final_to_best_error_ratio": final_to_best,
                "increase_count": reversals,
                "error_is_monotone": reversals == 0,
                "tail_relative_span": tail_span_relative,
            })
        summaries.append({
            "id": row["id"],
            "levels": row["refinement_levels"],
            "profiles": profile_summaries,
        })

    all_profiles = [profile for row in summaries for profile in row["profiles"]]
    checks = {
        "formulation_ids_unique": len({row["id"] for row in normalized}) == len(normalized),
        "refinement_levels_strictly_increase": all(
            all(a < b for a, b in zip(row["refinement_levels"], row["refinement_levels"][1:]))
            for row in normalized
        ),
        "force_errors_nonnegative": all(
            value >= 0.0
            for row in normalized
            for profile in row["error_profiles"]
            for value in profile
        ),
        "final_force_errors_meet_reference_band": all(
            profile["final_relative_error"] <= max_final_relative_error
            for profile in all_profiles
        ),
        "first_to_final_errors_improve": all(
            profile["initial_to_final_improvement"] >= min_initial_to_final_improvement
            for profile in all_profiles
        ),
        "final_errors_remain_near_best": all(
            profile["final_to_best_error_ratio"] <= max_final_to_best_error_ratio
            for profile in all_profiles
        ),
        "long_sweeps_reach_tail_plateau": all(
            profile["tail_relative_span"] is None
            or profile["tail_relative_span"] <= max_tail_relative_span
            for profile in all_profiles
        ),
    }
    return {
        "policy": "dual_formulation_force_error_convergence_gate_v1",
        "status": "ok" if all(checks.values()) else "needs_attention",
        "checks": checks,
        "issues": [name for name, ok in checks.items() if not ok],
        "reference_force": reference,
        "formulations": summaries,
        "metrics": {
            "formulation_count": len(summaries),
            "profile_count": len(all_profiles),
            "nonmonotone_profile_count": sum(
                not profile["error_is_monotone"] for profile in all_profiles
            ),
            "maximum_final_relative_error": max(
                profile["final_relative_error"] for profile in all_profiles
            ),
            "minimum_initial_to_final_improvement": min(
                profile["initial_to_final_improvement"] for profile in all_profiles
            ),
            "maximum_final_to_best_error_ratio": max(
                profile["final_to_best_error_ratio"] for profile in all_profiles
            ),
        },
        "lesson": (
            "Force-error convergence can oscillate. Gate the final analytic band, "
            "first-to-final improvement, final-to-best envelope, and tail plateau; "
            "record intermediate increases instead of rejecting them automatically."
        ),
    }
