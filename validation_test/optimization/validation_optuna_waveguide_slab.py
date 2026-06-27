"""Validation-class Optuna study for a waveguide dielectric slab.

The objective is a small CAE-style RF design task: choose a dielectric constant
and slab length in a WR-90-like TE10 guide so that the one-section slab is
nearly transparent at 10 GHz.  Analytic half-wave candidates are enqueued, then
Optuna records a reproducible TPE study.  The final feasible trial is selected
with radia-mcp's dependency-free record helpers.

Optuna is intentionally optional and external to radia-mcp.

Run:

    python validation_test/optimization/validation_optuna_waveguide_slab.py
"""

from __future__ import annotations

import datetime as dt
import importlib.metadata as metadata
import json
import math
import platform
import sys
from pathlib import Path

import optuna


HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
SRC = REPO / "packages" / "radia-mcp" / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from radia_mcp.radia_ngsolve.waveguide import waveguide_dielectric_slab_sparams  # noqa: E402
from radia_mcp.topology_optimization.global_optimizers import (  # noqa: E402
    best_feasible_record,
    constraint_violation,
)


OUT_JSON = HERE / "validation_optuna_waveguide_slab_summary.json"
FREQUENCY = 10.0e9
WIDTH_A = 0.02286
EPS_BOUNDS = (1.5, 4.0)
LENGTH_MM_BOUNDS = (2.0, 30.0)
DELIVERED_TARGET = 0.999
N_TRIALS = 40
SEED = 20260624


def _pkg_version(name: str) -> str | None:
    try:
        return metadata.version(name)
    except metadata.PackageNotFoundError:
        return None


def half_wave_length(eps_r: float) -> float:
    slab = waveguide_dielectric_slab_sparams(FREQUENCY, WIDTH_A, eps_r, 0.0)
    return math.pi / slab["beta_slab"]


def evaluate(eps_r: float, length_m: float) -> dict:
    slab = waveguide_dielectric_slab_sparams(FREQUENCY, WIDTH_A, eps_r, length_m)
    delivered = slab["S21_mag"] ** 2
    constraints = [DELIVERED_TARGET - delivered]
    return {
        "objective": slab["S11_mag"],
        "S11_mag": slab["S11_mag"],
        "S21_mag": slab["S21_mag"],
        "delivered_power_fraction": delivered,
        "unitarity": slab["unitarity"],
        "constraints": constraints,
        "constraint_violation": constraint_violation(constraints),
    }


def main() -> int:
    optuna.logging.set_verbosity(optuna.logging.WARNING)
    sampler = optuna.samplers.TPESampler(seed=SEED)
    study = optuna.create_study(direction="minimize", sampler=sampler)

    enqueued = []
    for eps_r in (1.5, 2.2, 3.0, 4.0):
        length_mm = 1.0e3 * half_wave_length(eps_r)
        if LENGTH_MM_BOUNDS[0] <= length_mm <= LENGTH_MM_BOUNDS[1]:
            params = {"eps_r": eps_r, "length_mm": length_mm}
            study.enqueue_trial(params)
            enqueued.append(params)

    def objective(trial):
        eps_r = trial.suggest_float("eps_r", *EPS_BOUNDS)
        length_mm = trial.suggest_float("length_mm", *LENGTH_MM_BOUNDS)
        result = evaluate(eps_r, length_mm * 1.0e-3)
        trial.set_user_attr("S21_mag", result["S21_mag"])
        trial.set_user_attr("delivered_power_fraction", result["delivered_power_fraction"])
        trial.set_user_attr("constraint_delivered", result["constraints"][0])
        trial.set_user_attr("constraint_violation", result["constraint_violation"])
        return result["objective"]

    study.optimize(objective, n_trials=N_TRIALS)

    records = []
    for trial in study.trials:
        constraints = [trial.user_attrs["constraint_delivered"]]
        records.append({
            "number": trial.number,
            "value": trial.value,
            "params": dict(trial.params),
            "constraints": constraints,
            "constraint_violation": constraint_violation(constraints),
            "delivered_power_fraction": trial.user_attrs["delivered_power_fraction"],
            "S21_mag": trial.user_attrs["S21_mag"],
        })
    best_feasible = best_feasible_record(records)
    best_trial = study.best_trial

    summary = {
        "schema": "radia.validation.optimization.optuna_waveguide_slab.v1",
        "generated_at_utc": dt.datetime.now(
            dt.timezone.utc
        ).isoformat(timespec="seconds").replace("+00:00", "Z"),
        "kind": "optuna_waveguide_slab_validation",
        "validation_class": True,
        "optuna_version": optuna.__version__,
        "versions": {
            "python_version": sys.version.split()[0],
            "platform": platform.platform(),
            "optuna_version": optuna.__version__,
            "radia_mcp_version": _pkg_version("radia-mcp"),
        },
        "seed": SEED,
        "n_trials": N_TRIALS,
        "frequency": FREQUENCY,
        "width_a": WIDTH_A,
        "eps_bounds": EPS_BOUNDS,
        "length_mm_bounds": LENGTH_MM_BOUNDS,
        "delivered_target": DELIVERED_TARGET,
        "enqueued": enqueued,
        "study_best": {
            "number": best_trial.number,
            "value": best_trial.value,
            "params": dict(best_trial.params),
        },
        "best_feasible": best_feasible,
        "records": records,
    }
    OUT_JSON.write_text(json.dumps(summary, indent=2), encoding="utf-8")

    print(f"[optuna] version={optuna.__version__}, trials={len(study.trials)}, seed={SEED}")
    print(f"[enqueue] {len(enqueued)} analytic half-wave candidates")
    print(
        f"[best] trial={best_trial.number} |S11|={best_trial.value:.6e} "
        f"eps={best_trial.params['eps_r']:.6g} "
        f"length={best_trial.params['length_mm']:.9f} mm"
    )
    print(
        f"[best feasible] trial={best_feasible['number']} "
        f"|S11|={best_feasible['value']:.6e} "
        f"delivered={best_feasible['delivered_power_fraction']:.12f} "
        f"violation={best_feasible['constraint_violation']:.3e}"
    )
    print(f"[OK] wrote {OUT_JSON}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
