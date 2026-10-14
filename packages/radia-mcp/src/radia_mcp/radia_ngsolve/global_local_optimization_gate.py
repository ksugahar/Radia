"""Solver-neutral global-search to local-polish replay validation."""
from __future__ import annotations
import json, math


def global_local_optimization_replay_gate(summary_json: str) -> dict[str, object]:
    try:
        row = json.loads(summary_json)
    except json.JSONDecodeError as exc:
        raise ValueError(f"summary_json must be valid JSON: {exc.msg}") from exc
    short = row.get("short_runs"); long = row.get("long_runs")
    if not isinstance(short, list) or not isinstance(long, list) or len(short) < 3 or len(short) != len(long):
        raise ValueError("short_runs and long_runs must contain the same three or more seeds")
    def finite(name, value):
        value=float(value)
        if not math.isfinite(value): raise ValueError(f"{name} must be finite")
        return value
    short_values=[finite("short best", r["best_f"]) for r in short]
    long_values=[finite("long best", r["best_f"]) for r in long]
    short_seeds=[int(r["seed"]) for r in short]; long_seeds=[int(r["seed"]) for r in long]
    analytic=finite("analytic minimum", row["analytic_minimum"])
    global_best=finite("independent global best", row["independent_global_best_f"])
    polished=finite("polished best", row["polished_best_f"])
    checks={
        "source_objective_replayed": finite("source objective error", row["source_objective_max_abs_error"]) <= 1e-12,
        "same_unique_seed_set": short_seeds == long_seeds and len(set(long_seeds)) == len(long_seeds),
        "short_histories_monotone": all(r.get("history_monotone") is True for r in short),
        "long_histories_monotone": all(r.get("history_monotone") is True for r in long),
        "long_budget_not_worse": all(b <= a + 1e-12 for a,b in zip(short_values,long_values)),
        "multi_seed_global_basin_resolved": min(long_values) <= analytic + 1e-4,
        "independent_global_solver_agrees": global_best <= analytic + 1e-4,
        "local_polish_matches_analytic_minimum": abs(polished-analytic) <= 1e-12,
        "polished_gradient_small": finite("gradient norm", row["polished_gradient_norm"]) <= 1e-6,
        "central_gradient_matches": finite("central gradient error", row["central_gradient_relative_error"]) <= 1e-7,
        "complex_step_gradient_matches": finite("complex-step error", row["complex_step_gradient_relative_error"]) <= 1e-12,
    }
    return {"policy":"global_local_optimization_replay_gate_v1","status":"ok" if all(checks.values()) else "needs_attention",
            "checks":checks,"issues":[k for k,v in checks.items() if not v],
            "metrics":{"seed_count":len(long),"short_best_f":min(short_values),"long_best_f":min(long_values),
                       "independent_global_best_f":global_best,"polished_best_f":polished},
            "lesson":"Treat stochastic global search as a budgeted basin finder. Require seed replay, monotone histories, an independent global solver, local polish, and derivative checks before accepting the optimum."}
