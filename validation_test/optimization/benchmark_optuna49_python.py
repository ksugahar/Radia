"""Benchmark the pinned upstream Optuna 4.9.0 TPE oracle.

Run this and ``benchmark_matlab_optuna49.m`` on the same otherwise-idle host.
The first three repeats are warm-ups; the reported value is the median of the
remaining eight repeats. The workloads intentionally exclude persistence and
parallel scheduling because those are MATLAB extensions rather than shared
Optuna behavior.
"""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
import json
import os
from pathlib import Path
import platform
import statistics
import sys
import time
import warnings

import numpy
import optuna
import scipy


TRIALS = 100
REPEATS = 11
WARMUP_REPEATS = 3
DATAFRAME_TRIALS = 1000
EXPECTED_SCALAR_CHECKSUM = 20.040135043951892
EXPECTED_GROUPED_CHECKSUM = 104.33176385944043


def _scalar() -> float:
    study = optuna.create_study(
        sampler=optuna.samplers.TPESampler(seed=37, n_startup_trials=4)
    )
    checksum = 0.0
    for _ in range(TRIALS):
        trial = study.ask()
        x = trial.suggest_float("x", -2.0, 2.0)
        study.tell(trial, (x - 0.25) ** 2)
        checksum += x
    return checksum


def _grouped() -> float:
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        sampler = optuna.samplers.TPESampler(
            seed=101,
            n_startup_trials=4,
            multivariate=True,
            group=True,
        )
    study = optuna.create_study(sampler=sampler)
    checksum = 0.0
    for _ in range(TRIALS):
        trial = study.ask()
        branch = trial.suggest_categorical("branch", ["left", "right"])
        x = trial.suggest_float("x", -1.0, 1.0)
        if branch == "left":
            y = trial.suggest_float("y", 0.0, 2.0)
            value = (x - 0.2) ** 2 + (y - 0.4) ** 2
            checksum += x + y
        else:
            z = trial.suggest_int("z", 1, 5)
            value = (x + 0.1) ** 2 + 0.05 * z
            checksum += x + z
        study.tell(trial, value)
    return checksum


def _measure(workload, expected_checksum: float) -> dict[str, object]:
    durations: list[float] = []
    checksums: list[float] = []
    for _ in range(REPEATS):
        started = time.perf_counter()
        checksums.append(workload())
        durations.append(time.perf_counter() - started)
    if any(abs(value - expected_checksum) > 1e-12 for value in checksums):
        raise RuntimeError(
            f"Seeded Optuna checksum changed: {checksums[-1]!r} != "
            f"{expected_checksum!r}"
        )
    median_seconds = statistics.median(durations[WARMUP_REPEATS:])
    return {
        "all_seconds": durations,
        "median_warmed_seconds": median_seconds,
        "trials_per_second": TRIALS / median_seconds,
        "checksum": checksums[-1],
    }


def _measure_trials_dataframe() -> dict[str, object]:
    distributions = {
        "x": optuna.distributions.FloatDistribution(0.0, 1.0),
        "mode": optuna.distributions.CategoricalDistribution(["A", "B"]),
    }
    study = optuna.create_study()
    for index in range(DATAFRAME_TRIALS):
        study.add_trial(
            optuna.trial.create_trial(
                value=float(index),
                params={
                    "x": index / DATAFRAME_TRIALS,
                    "mode": ["A", "B"][index % 2],
                },
                distributions=distributions,
                user_attrs={"owner": "lab"},
                system_attrs={"origin": "benchmark"},
            )
        )
    attrs = (
        "number",
        "value",
        "params",
        "user_attrs",
        "system_attrs",
        "state",
    )
    durations: list[float] = []
    frame = None
    for _ in range(REPEATS):
        started = time.perf_counter()
        frame = study.trials_dataframe(attrs=attrs)
        durations.append(time.perf_counter() - started)
    assert frame is not None
    expected_columns = [
        "number",
        "value",
        "params_mode",
        "params_x",
        "user_attrs_owner",
        "system_attrs_origin",
        "state",
    ]
    if list(frame.columns) != expected_columns or len(frame) != DATAFRAME_TRIALS:
        raise RuntimeError("trials_dataframe benchmark contract changed")
    median_seconds = statistics.median(durations[WARMUP_REPEATS:])
    return {
        "all_seconds": durations,
        "median_warmed_seconds": median_seconds,
        "rows_per_second": DATAFRAME_TRIALS / median_seconds,
        "rows": DATAFRAME_TRIALS,
        "columns": len(expected_columns),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    if optuna.__version__ != "4.9.0":
        raise RuntimeError(
            f"This benchmark requires optuna==4.9.0, found {optuna.__version__}"
        )
    optuna.logging.set_verbosity(optuna.logging.ERROR)
    result = {
        "schema": "radia.validation.optuna49-performance-runtime.v1",
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "runtime": "python-upstream",
        "host": os.environ.get("COMPUTERNAME", platform.node()),
        "versions": {
            "python": sys.version.split()[0],
            "optuna": optuna.__version__,
            "numpy": numpy.__version__,
            "scipy": scipy.__version__,
        },
        "settings": {
            "trials": TRIALS,
            "total_repeats": REPEATS,
            "warmup_repeats": WARMUP_REPEATS,
            "reported_repeats": REPEATS - WARMUP_REPEATS,
        },
        "scalar": _measure(_scalar, EXPECTED_SCALAR_CHECKSUM),
        "grouped_conditional": _measure(
            _grouped, EXPECTED_GROUPED_CHECKSUM
        ),
        "trials_dataframe": _measure_trials_dataframe(),
    }
    encoded = json.dumps(result, indent=2) + "\n"
    if args.output:
        args.output.write_text(encoded, encoding="utf-8")
    print(encoded, end="")


if __name__ == "__main__":
    main()
