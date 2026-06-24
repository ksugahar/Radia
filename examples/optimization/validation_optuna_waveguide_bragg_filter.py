"""Validation-class Optuna study for a waveguide Bragg stopband filter.

The objective is a compact CAE-style RF optimization task: tune a 3-period
high/low-permittivity TE10 waveguide stack to suppress transmission at 10 GHz.
Analytic quarter-wave candidates are enqueued, then Optuna records a
reproducible TPE study.  The final feasible trial is selected with radia-mcp's
dependency-free record helpers.

Optuna is intentionally optional and external to radia-mcp.

Run:

    python examples/optimization/validation_optuna_waveguide_bragg_filter.py
"""

from __future__ import annotations

import json
import math
import sys
from pathlib import Path

import optuna


HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
SRC = REPO / "packages" / "radia-mcp" / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from radia_mcp.radia_ngsolve.waveguide import C0, waveguide_cascade_sparams  # noqa: E402
from radia_mcp.topology_optimization.global_optimizers import (  # noqa: E402
    best_feasible_record,
    constraint_violation,
)


OUT_JSON = HERE / "validation_optuna_waveguide_bragg_filter_summary.json"
FREQUENCY = 10.0e9
WIDTH_A = 0.02286
N_PERIODS = 3
REFLECTED_TARGET = 0.95
N_TRIALS = 50
SEED = 20260624

BOUNDS = {
    "eps_hi": (3.0, 6.0),
    "eps_lo": (1.05, 1.8),
    "scale_hi": (0.85, 1.15),
    "scale_lo": (0.85, 1.15),
}

ENQUEUED = [
    {"eps_hi": 4.0, "eps_lo": 1.2, "scale_hi": 1.0, "scale_lo": 1.0},
    {"eps_hi": 6.0, "eps_lo": 1.05, "scale_hi": 1.0, "scale_lo": 1.0},
    {"eps_hi": 5.0, "eps_lo": 1.1, "scale_hi": 1.0, "scale_lo": 1.0},
]


def quarter_wave_length(eps_r: float) -> float:
    k0 = 2.0 * math.pi * FREQUENCY / C0
    kc = math.pi / WIDTH_A
    return math.pi / (2.0 * math.sqrt(eps_r * k0 * k0 - kc * kc))


def sections_from_params(params: dict) -> list[tuple[float, float]]:
    hi = quarter_wave_length(params["eps_hi"]) * params["scale_hi"]
    lo = quarter_wave_length(params["eps_lo"]) * params["scale_lo"]
    sections = []
    for _ in range(N_PERIODS):
        sections.append((hi, params["eps_hi"]))
        sections.append((lo, params["eps_lo"]))
    return sections


def evaluate(params: dict) -> dict:
    sections = sections_from_params(params)
    s = waveguide_cascade_sparams(FREQUENCY, WIDTH_A, sections)
    reflected = s["S11_mag"] ** 2
    constraints = [REFLECTED_TARGET - reflected]
    return {
        "objective": s["S21_mag"],
        "S11_mag": s["S11_mag"],
        "S21_mag": s["S21_mag"],
        "reflected_power": reflected,
        "transmitted_power": s["S21_mag"] ** 2,
        "unitarity": s["unitarity"],
        "max_section_length_mm": 1.0e3 * max(length for length, _ in sections),
        "sections_mm": [1.0e3 * length for length, _ in sections],
        "constraints": constraints,
        "constraint_violation": constraint_violation(constraints),
    }


def _trial_params(trial) -> dict:
    params = {
        "eps_hi": trial.suggest_float("eps_hi", *BOUNDS["eps_hi"]),
        "eps_lo": trial.suggest_float("eps_lo", *BOUNDS["eps_lo"]),
        "scale_hi": trial.suggest_float("scale_hi", *BOUNDS["scale_hi"]),
        "scale_lo": trial.suggest_float("scale_lo", *BOUNDS["scale_lo"]),
    }
    if params["eps_hi"] <= params["eps_lo"]:
        # Outside the intended high/low stack ordering; keep it unattractive.
        params["eps_hi"], params["eps_lo"] = params["eps_lo"], params["eps_hi"]
    return params


def main() -> int:
    optuna.logging.set_verbosity(optuna.logging.WARNING)
    sampler = optuna.samplers.TPESampler(seed=SEED)
    study = optuna.create_study(direction="minimize", sampler=sampler)

    for params in ENQUEUED:
        study.enqueue_trial(params)

    def objective(trial):
        params = _trial_params(trial)
        result = evaluate(params)
        for key in ("S11_mag", "S21_mag", "reflected_power", "transmitted_power",
                    "unitarity", "max_section_length_mm", "constraint_violation"):
            trial.set_user_attr(key, result[key])
        trial.set_user_attr("constraint_reflected", result["constraints"][0])
        trial.set_user_attr("sections_mm", result["sections_mm"])
        return result["objective"]

    study.optimize(objective, n_trials=N_TRIALS)

    records = []
    for trial in study.trials:
        constraints = [trial.user_attrs["constraint_reflected"]]
        records.append({
            "number": trial.number,
            "value": trial.value,
            "params": dict(trial.params),
            "constraints": constraints,
            "constraint_violation": constraint_violation(constraints),
            "S11_mag": trial.user_attrs["S11_mag"],
            "S21_mag": trial.user_attrs["S21_mag"],
            "reflected_power": trial.user_attrs["reflected_power"],
            "transmitted_power": trial.user_attrs["transmitted_power"],
            "unitarity": trial.user_attrs["unitarity"],
            "max_section_length_mm": trial.user_attrs["max_section_length_mm"],
            "sections_mm": trial.user_attrs["sections_mm"],
        })

    best_feasible = best_feasible_record(records)
    best_trial = study.best_trial
    enqueued_records = [records[i] for i in range(len(ENQUEUED))]
    max_unitarity_error = max(abs(r["unitarity"] - 1.0) for r in records)

    summary = {
        "kind": "optuna_waveguide_bragg_filter_validation",
        "validation_class": True,
        "optuna_version": optuna.__version__,
        "seed": SEED,
        "n_trials": N_TRIALS,
        "frequency": FREQUENCY,
        "width_a": WIDTH_A,
        "n_periods": N_PERIODS,
        "bounds": BOUNDS,
        "reflected_target": REFLECTED_TARGET,
        "enqueued": ENQUEUED,
        "enqueued_records": enqueued_records,
        "study_best": {
            "number": best_trial.number,
            "value": best_trial.value,
            "params": dict(best_trial.params),
        },
        "best_feasible": best_feasible,
        "max_unitarity_error": max_unitarity_error,
        "records": records,
    }

    assert enqueued_records[0]["constraint_violation"] == 0.0
    assert best_feasible["constraint_violation"] == 0.0
    assert best_feasible["reflected_power"] >= REFLECTED_TARGET
    assert best_feasible["S21_mag"] < 0.10
    assert best_feasible["value"] <= enqueued_records[0]["value"]
    assert max_unitarity_error < 1.0e-12

    OUT_JSON.write_text(json.dumps(summary, indent=2), encoding="utf-8")

    print(f"[optuna] version={optuna.__version__}, trials={len(study.trials)}, seed={SEED}")
    print(f"[enqueue] {len(ENQUEUED)} analytic quarter-wave candidates")
    for rec in enqueued_records:
        print(
            f"  trial={rec['number']:2d} |S21|={rec['S21_mag']:.6f} "
            f"reflected={rec['reflected_power']:.6f} "
            f"eps_hi={rec['params']['eps_hi']:.4g} eps_lo={rec['params']['eps_lo']:.4g}"
        )
    print(
        f"[best feasible] trial={best_feasible['number']} "
        f"|S21|={best_feasible['S21_mag']:.6f} "
        f"|S11|={best_feasible['S11_mag']:.6f} "
        f"reflected={best_feasible['reflected_power']:.12f} "
        f"violation={best_feasible['constraint_violation']:.3e}"
    )
    print(f"[OK] wrote {OUT_JSON}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
