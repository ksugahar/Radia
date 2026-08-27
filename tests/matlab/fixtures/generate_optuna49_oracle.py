"""Generate deterministic MATLAB parity fixtures from upstream Optuna 4.9.0."""

from __future__ import annotations

import importlib.metadata
import inspect
import json
import logging
import sys
import warnings
from pathlib import Path

import numpy as np
import optuna
import optuna.importance
import optuna.terminator
import scipy
from optuna.trial import TrialState

EXPECTED_VERSION = "4.9.0"


def _exception_contract() -> dict[str, object]:
    classes = [
        "CLIUsageError",
        "DuplicatedStudyError",
        "ExperimentalWarning",
        "OptunaError",
        "StorageInternalError",
        "TrialPruned",
        "UpdateFinishedTrialError",
    ]
    result: dict[str, object] = {}
    for name in classes:
        exception_type = getattr(optuna.exceptions, name)
        cases = []
        for arguments in [(), ("oracle message",), ("a", "b")]:
            exception = exception_type(*arguments)
            exception.add_note("oracle note")
            cases.append(
                {
                    "args": list(exception.args),
                    "message": str(exception),
                    "notes": list(exception.__notes__),
                    "with_traceback_identity": exception.with_traceback(None)
                    is exception,
                }
            )
        result[name] = {
            "cases": cases,
            "is_optuna_error": issubclass(exception_type, optuna.exceptions.OptunaError),
            "is_warning": issubclass(exception_type, Warning),
        }
    return result


def _logging_contract() -> dict[str, object]:
    module = optuna.logging
    constants = {
        name: getattr(module, name)
        for name in ["CRITICAL", "DEBUG", "ERROR", "FATAL", "INFO", "WARN", "WARNING"]
    }
    module.set_verbosity(module.INFO)
    initial = module.get_verbosity()
    logger = module.get_logger("unit")
    formatter = module.create_default_formatter()
    module.set_verbosity(module.DEBUG)
    after_set = module.get_verbosity()
    module.disable_default_handler()
    module.disable_propagation()
    root = module.get_logger("optuna")
    disabled = {"handlers": len(root.handlers), "propagate": root.propagate}
    module.enable_default_handler()
    module.enable_propagation()
    enabled = {"handlers": len(root.handlers), "propagate": root.propagate}
    module.set_verbosity(module.WARNING)
    module.disable_propagation()
    return {
        "after_set": after_set,
        "constants": constants,
        "disabled": disabled,
        "enabled": enabled,
        "formatter": {"date_format": formatter.datefmt, "format": formatter._style._fmt},
        "initial": initial,
        "logger": {
            "handlers": len(logger.handlers),
            "level": logger.level,
            "name": logger.name,
            "propagate": logger.propagate,
        },
    }


def _sampler_seed_default_contract() -> dict[str, object]:
    constructors = [
        "RandomSampler",
        "TPESampler",
        "CmaEsSampler",
        "GPSampler",
        "GridSampler",
        "NSGAIISampler",
        "NSGAIIISampler",
        "QMCSampler",
        "BruteForceSampler",
    ]
    defaults: dict[str, object] = {}
    for name in constructors:
        parameter = inspect.signature(getattr(optuna.samplers, name)).parameters["seed"]
        if parameter.default is not None:
            raise RuntimeError(f"Optuna 4.9.0 {name}.seed no longer defaults to None.")
        defaults[name] = {
            "parameter": parameter.name,
            "default_is_none": True,
        }
    return {
        "constructors": defaults,
        "semantic": "fresh nondeterministic entropy per unseeded sampler instance",
        "exact_sequence_oracle_requires_explicit_seed": True,
    }


def _sampler_public_member_contract() -> dict[str, object]:
    def make_study() -> tuple[optuna.Study, optuna.trial.FrozenTrial, optuna.trial.FrozenTrial]:
        study = optuna.create_study()
        for value in [0.2, 0.1]:
            trial = study.ask()
            trial.suggest_float("x", 0.0, 1.0, step=0.1)
            trial.suggest_int("y", 1, 5)
            study.tell(trial, value)
        study.ask()
        trials = study.get_trials(deepcopy=False)
        return study, trials[-1], trials[0]

    samplers = {
        "BruteForceSampler": optuna.samplers.BruteForceSampler(seed=7),
        "CmaEsSampler": optuna.samplers.CmaEsSampler(
            seed=7, n_startup_trials=99
        ),
        "GPSampler": optuna.samplers.GPSampler(seed=7, n_startup_trials=99),
        "GridSampler": optuna.samplers.GridSampler(
            {"x": [0.0, 1.0], "y": [1, 3]}, seed=7
        ),
        "NSGAIIISampler": optuna.samplers.NSGAIIISampler(
            seed=7, population_size=4
        ),
        "NSGAIISampler": optuna.samplers.NSGAIISampler(
            seed=7, population_size=4
        ),
        "PartialFixedSampler": optuna.samplers.PartialFixedSampler(
            {"x": 0.25},
            optuna.samplers.TPESampler(
                seed=7, n_startup_trials=1, multivariate=True
            ),
        ),
        "QMCSampler": optuna.samplers.QMCSampler(
            seed=7, warn_asynchronous_seeding=False
        ),
        "RandomSampler": optuna.samplers.RandomSampler(seed=7),
        "TPESampler": optuna.samplers.TPESampler(
            seed=7, n_startup_trials=1, multivariate=True
        ),
    }
    sampler_contract: dict[str, object] = {}
    for name, sampler in samplers.items():
        study, running, complete = make_study()
        search_space = sampler.infer_relative_search_space(study, running)
        before_result = sampler.before_trial(study, running)
        after_result = sampler.after_trial(
            study, complete, TrialState.COMPLETE, [complete.value]
        )
        sampler_contract[name] = {
            "after_trial_returns_none": after_result is None,
            "before_trial_returns_none": before_result is None,
            "infer_relative_search_space_keys": sorted(search_space),
        }

    grid = optuna.samplers.GridSampler({"x": [0, 1]}, seed=7)
    grid_study = optuna.create_study(sampler=grid)
    grid_exhausted_before = grid.is_exhausted(grid_study)
    grid_study.optimize(lambda trial: float(trial.suggest_int("x", 0, 1)), n_trials=2)

    relative_study = optuna.create_study(
        sampler=optuna.samplers.TPESampler(
            seed=7, n_startup_trials=0, multivariate=True
        )
    )
    complete_relative = relative_study.ask()
    complete_relative.suggest_float("x", 0.0, 1.0, step=0.1)
    complete_relative.suggest_int("y", 1, 5)
    relative_study.tell(complete_relative, 0.2)
    running_relative = relative_study.ask()

    crossovers = {
        "blxalpha": optuna.samplers.nsgaii.BLXAlphaCrossover(),
        "sbx": optuna.samplers.nsgaii.SBXCrossover(),
        "spx": optuna.samplers.nsgaii.SPXCrossover(),
        "undx": optuna.samplers.nsgaii.UNDXCrossover(),
        "uniform": optuna.samplers.nsgaii.UniformCrossover(),
        "vsbx": optuna.samplers.nsgaii.VSBXCrossover(),
    }
    try:
        optuna.samplers.nsgaii.BaseCrossover()
    except Exception as error:  # noqa: BLE001 - public abstract-class contract.
        base_crossover_error = type(error).__name__
    else:
        raise RuntimeError("Optuna BaseCrossover unexpectedly became concrete.")

    mapped_namespaces: dict[str, str] = {}
    for name in [
        "distributions",
        "exceptions",
        "importance",
        "pruners",
        "samplers",
        "search_space",
        "storages",
        "study",
        "trial",
    ]:
        module = getattr(optuna, name)
        if not inspect.ismodule(module):
            raise RuntimeError(f"Optuna public namespace {name!r} is not a module.")
        mapped_namespaces[name] = module.__name__

    fixed = optuna.trial.FixedTrial({"x": 0.5})
    return {
        "base_crossover_instantiation_error": base_crossover_error,
        "crossover_n_parents": {
            name: crossover.n_parents for name, crossover in crossovers.items()
        },
        "fixed_trial_datetime_start_is_not_none": fixed.datetime_start is not None,
        "grid_is_exhausted_after": grid.is_exhausted(grid_study),
        "grid_is_exhausted_before": grid_exhausted_before,
        "mapped_namespaces": mapped_namespaces,
        "nsgaii_population_size": optuna.samplers.NSGAIISampler(
            population_size=5
        ).population_size,
        "nsgaiii_population_size": optuna.samplers.NSGAIIISampler(
            population_size=6
        ).population_size,
        "samplers": sampler_contract,
        "trial_datetime_start_is_not_none": running_relative.datetime_start is not None,
        "trial_relative_params_keys": sorted(running_relative.relative_params),
    }


def _random_trials() -> list[dict[str, object]]:
    study = optuna.create_study(sampler=optuna.samplers.RandomSampler(seed=123))
    rows: list[dict[str, object]] = []
    for _ in range(16):
        trial = study.ask()
        row = {
            "x": trial.suggest_float("x", -1.0, 1.0),
            "q": trial.suggest_float("q", 0.0, 1.0, step=0.2),
            "mesh": trial.suggest_int("mesh", 1, 9, step=2),
            "log_mesh": trial.suggest_int("log_mesh", 1, 100, log=True),
            "mode": trial.suggest_categorical("mode", ["A", "B", "C"]),
        }
        study.tell(trial, row["x"])
        rows.append(row)
    return rows


def _tpe_trials() -> list[float]:
    study = optuna.create_study(
        sampler=optuna.samplers.TPESampler(seed=37, n_startup_trials=4)
    )
    values: list[float] = []
    for _ in range(40):
        trial = study.ask()
        x = trial.suggest_float("x", -2.0, 2.0)
        study.tell(trial, (x - 0.25) ** 2)
        values.append(x)
    return values


def _custom_tpe_trials() -> dict[str, object]:
    def gamma(count: int) -> int:
        return min(3, count)

    def weights(count: int) -> np.ndarray:
        return np.linspace(0.2, 1.0, count) if count else np.asarray([])

    single = optuna.create_study(
        sampler=optuna.samplers.TPESampler(
            seed=97,
            n_startup_trials=4,
            gamma=gamma,
            weights=weights,
        )
    )
    single_values: list[float] = []
    for _ in range(32):
        trial = single.ask()
        x = trial.suggest_float("x", -2.0, 2.0)
        single.tell(trial, (x - 0.35) ** 2)
        single_values.append(x)

    multi = optuna.create_study(
        directions=["minimize", "minimize"],
        sampler=optuna.samplers.TPESampler(
            seed=101,
            n_startup_trials=4,
            gamma=gamma,
            weights=weights,
        ),
    )
    multi_values: list[dict[str, float]] = []
    for _ in range(24):
        trial = multi.ask()
        x = trial.suggest_float("x", -2.0, 2.0)
        y = trial.suggest_float("y", -1.0, 3.0)
        multi.tell(trial, [(x - 0.4) ** 2 + 0.1 * y * y, (y + 0.2) ** 2])
        multi_values.append({"x": x, "y": y})

    return {"single": single_values, "multi": multi_values}


def _tpe_group_contract() -> dict[str, object]:
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        sampler = optuna.samplers.TPESampler(
            seed=101,
            n_startup_trials=4,
            multivariate=True,
            group=True,
        )
    study = optuna.create_study(sampler=sampler)
    rows: list[dict[str, object]] = []
    for _ in range(12):
        trial = study.ask()
        branch = trial.suggest_categorical("branch", ["left", "right"])
        x = trial.suggest_float("x", -1.0, 1.0)
        row: dict[str, object] = {
            "number": trial.number,
            "branch": branch,
            "x": x,
            "y": None,
            "z": None,
        }
        if branch == "left":
            y = trial.suggest_float("y", 0.0, 2.0)
            row["y"] = y
            value = (x - 0.2) ** 2 + (y - 0.4) ** 2
        else:
            z = trial.suggest_int("z", 1, 5)
            row["z"] = z
            value = (x + 0.1) ** 2 + 0.05 * z
        study.tell(trial, value)
        rows.append(row)

    logger = logging.getLogger("optuna.samplers._tpe.sampler")

    class Records(logging.Handler):
        def __init__(self) -> None:
            super().__init__()
            self.messages: list[str] = []

        def emit(self, record: logging.LogRecord) -> None:
            self.messages.append(record.getMessage())

    def capture(*, group: bool, warn: bool) -> list[str]:
        handler = Records()
        logger.addHandler(handler)
        try:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                local_sampler = optuna.samplers.TPESampler(
                    seed=103,
                    n_startup_trials=1,
                    multivariate=True,
                    group=group,
                    warn_independent_sampling=warn,
                )
            local_study = optuna.create_study(sampler=local_sampler)
            first = local_study.ask()
            x = first.suggest_float("x", 0.0, 1.0)
            local_study.tell(first, x)
            second = local_study.ask()
            x = second.suggest_float("x", 0.0, 1.0)
            y = second.suggest_float("y", 0.0, 1.0)
            local_study.tell(second, x + y)
            third = local_study.ask()
            y = third.suggest_float("y", 0.0, 1.0)
            local_study.tell(third, y)
        finally:
            logger.removeHandler(handler)
        return handler.messages

    return {
        "sequence": rows,
        "independent_warning_enabled_count": len(
            capture(group=False, warn=True)
        ),
        "independent_warning_disabled_count": len(
            capture(group=False, warn=False)
        ),
        "group_warning_count": len(capture(group=True, warn=True)),
    }


def _tpe_categorical_distance_trials() -> list[dict[str, object]]:
    levels = ["zero", "one", "two", "three"]
    positions = {level: index for index, level in enumerate(levels)}

    def distance(first: str, second: str) -> float:
        return float(abs(positions[first] - positions[second]))

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        sampler = optuna.samplers.TPESampler(
            seed=107,
            n_startup_trials=4,
            categorical_distance_func={"level": distance},
        )
    study = optuna.create_study(sampler=sampler)
    rows: list[dict[str, object]] = []
    for _ in range(18):
        trial = study.ask()
        level = trial.suggest_categorical("level", levels)
        position = positions[level]
        study.tell(trial, (position - 1.3) ** 2)
        rows.append(
            {"number": trial.number, "level": level, "position": position}
        )
    return rows


def _multiobjective_tpe_trials() -> list[dict[str, float]]:
    study = optuna.create_study(
        directions=["minimize", "minimize"],
        sampler=optuna.samplers.TPESampler(seed=41, n_startup_trials=4),
    )
    rows: list[dict[str, float]] = []
    for _ in range(32):
        trial = study.ask()
        x = trial.suggest_float("x", -2.0, 2.0)
        y = trial.suggest_float("y", -1.0, 3.0)
        study.tell(trial, [(x - 0.4) ** 2 + 0.1 * y * y, (y + 0.2) ** 2 + 0.1 * x * x])
        rows.append({"x": x, "y": y})
    return rows


def _mixed_tpe_trials() -> list[dict[str, object]]:
    study = optuna.create_study(
        sampler=optuna.samplers.TPESampler(seed=43, n_startup_trials=4)
    )
    rows: list[dict[str, object]] = []
    # At trial 15 this deliberately small categorical history reaches an
    # exact likelihood-ratio tie. NumPy/libm then select different categories
    # by one ULP on different CPU backends, which is not a stable Optuna
    # contract. Keep this mixed-distribution oracle before that backend tie;
    # longer scalar and multivariate TPE sequences are covered separately.
    for _ in range(14):
        trial = study.ask()
        row = {
            "x": trial.suggest_float("x", -1.0, 1.0),
            "q": trial.suggest_float("q", 0.0, 1.0, step=0.2),
            "mesh": trial.suggest_int("mesh", 1, 9, step=2),
            "log_mesh": trial.suggest_int("log_mesh", 1, 100, log=True),
            "mode": trial.suggest_categorical("mode", ["A", "B", "C"]),
        }
        loss = (
            (float(row["x"]) - 0.2) ** 2
            + (float(row["q"]) - 0.6) ** 2
            + 0.01 * float(row["mesh"])
            + 0.001 * float(row["log_mesh"])
            + {"A": 0.2, "B": 0.0, "C": 0.35}[str(row["mode"])]
        )
        study.tell(trial, loss)
        rows.append(row)
    return rows


def _random_state_contract() -> dict[str, object]:
    uniform_rng = np.random.RandomState(37)
    uniforms = uniform_rng.rand(400)
    normal_rng = np.random.RandomState(37)
    permutation_rng = np.random.RandomState(123)
    integer_rng = np.random.RandomState(123)
    positions = [0, 1, 2, 112, 113, 119, 120, 311, 312, 399]
    return {
        "uniform_positions_zero_based": positions,
        "uniform_values": [float(uniforms[index]) for index in positions],
        "normal_values": normal_rng.randn(7).tolist(),
        "permutation_one_based": (permutation_rng.permutation(10) + 1).tolist(),
        "integers_one_based": integer_rng.randint(1, 11, size=8).tolist(),
    }


def _core_api_contract() -> dict[str, object]:
    queued = optuna.create_study(sampler=optuna.samplers.RandomSampler(seed=7))
    params = {"x": 0.25, "mesh": 5, "mode": "B"}
    queued.enqueue_trial(params, user_attrs={"source": "baseline"})
    queued.enqueue_trial(params, skip_if_exists=True)
    waiting_count = len(queued.trials)
    queued_trial = queued.ask()
    queued_values = {
        "number": queued_trial.number,
        "x": queued_trial.suggest_float("x", 0.0, 1.0),
        "mesh": queued_trial.suggest_int("mesh", 1, 9, step=2),
        "mode": queued_trial.suggest_categorical("mode", ["A", "B"]),
        "source": queued_trial.user_attrs["source"],
    }
    queued.tell(queued_trial, 1.0)

    invalid = optuna.create_study(sampler=optuna.samplers.RandomSampler(seed=9))
    invalid.enqueue_trial({"x": 2.0})
    invalid_trial = invalid.ask()
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        invalid_fallback = invalid_trial.suggest_float("x", 0.0, 1.0)

    fixed = optuna.create_study(sampler=optuna.samplers.RandomSampler(seed=19))
    fixed_trial = fixed.ask(
        {
            "x": optuna.distributions.FloatDistribution(-1.0, 1.0),
            "mesh": optuna.distributions.IntDistribution(1, 9, step=2),
            "mode": optuna.distributions.CategoricalDistribution(["A", "B"]),
        }
    )

    return {
        "enqueue": {
            "waiting_count_after_skip": waiting_count,
            "values": queued_values,
            "final_state": queued.trials[0].state.name,
        },
        "invalid_enqueued_fallback": float(invalid_fallback),
        "ask_fixed_distributions": {
            "x": float(fixed_trial.params["x"]),
            "mesh": int(fixed_trial.params["mesh"]),
            "mode": fixed_trial.params["mode"],
        },
    }


def _distribution_contract() -> dict[str, object]:
    study = optuna.create_study(sampler=optuna.samplers.RandomSampler(seed=4))
    trial = study.ask()
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        first = trial.suggest_int("mesh order", 2, 10, step=3)
        repeated = trial.suggest_int("mesh order", 2, 10, step=2)
    errors: dict[str, str] = {}
    for name, call in {
        "different_kind": lambda: trial.suggest_float("mesh order", 2.0, 10.0),
        "different_log": lambda: trial.suggest_int("mesh order", 2, 10, log=True),
    }.items():
        try:
            call()
        except Exception as error:  # noqa: BLE001 - fixture records public error type.
            errors[name] = type(error).__name__

    categorical = study.ask()
    categorical.suggest_categorical("choice", ["A", "B"])
    try:
        categorical.suggest_categorical("choice", ["A", "C"])
    except Exception as error:  # noqa: BLE001
        errors["dynamic_categorical"] = type(error).__name__

    collision = study.ask()
    first_name = collision.suggest_float("a-b", 0.0, 1.0)
    second_name = collision.suggest_float("a_b", 0.0, 1.0)

    stepped = study.ask()
    with warnings.catch_warnings(record=True) as stepped_warnings:
        warnings.simplefilter("always")
        stepped_value = stepped.suggest_float("q", 0.0, 1.0, step=0.3)
    distribution = stepped.distributions["q"]
    return {
        "integer_value": int(first),
        "inconsistent_repeat_value": int(repeated),
        "integer_warning_categories": [type(item.message).__name__ for item in caught],
        "errors": errors,
        "colliding_names": {"a-b": float(first_name), "a_b": float(second_name)},
        "stepped_float": {
            "value": float(stepped_value),
            "effective_high": float(distribution.high),
            "warning_categories": [
                type(item.message).__name__ for item in stepped_warnings
            ],
        },
    }


def _search_space_contract() -> dict[str, object]:
    distributions = optuna.distributions
    base = {
        "x": distributions.FloatDistribution(0.0, 1.0),
        "fixed": distributions.FloatDistribution(2.0, 2.0),
        "cat": distributions.CategoricalDistribution(["A", "B"]),
    }
    values = {"x": 0.5, "fixed": 2.0, "cat": "A", "z": 1}

    def frozen(state: TrialState, current: dict[str, object]) -> optuna.trial.FrozenTrial:
        return optuna.trial.create_trial(
            state=state,
            value=1.0 if state == TrialState.COMPLETE else None,
            params={name: values[name] for name in current},
            distributions=current,
        )

    trials = [
        frozen(TrialState.COMPLETE, base),
        frozen(
            TrialState.COMPLETE,
            {**base, "z": distributions.IntDistribution(1, 3)},
        ),
        frozen(
            TrialState.PRUNED,
            {
                "x": distributions.FloatDistribution(-1.0, 1.0),
                "fixed": base["fixed"],
                "cat": base["cat"],
            },
        ),
        frozen(
            TrialState.FAIL,
            {"x": distributions.FloatDistribution(-2.0, 2.0)},
        ),
        frozen(
            TrialState.WAITING,
            {"x": distributions.FloatDistribution(-3.0, 3.0)},
        ),
    ]
    without_pruned = optuna.search_space.intersection_search_space(trials)
    with_pruned = optuna.search_space.intersection_search_space(
        trials, include_pruned=True
    )

    study = optuna.create_study()
    study.add_trials(trials)
    calculator = optuna.search_space.IntersectionSearchSpace()
    calculated = calculator.calculate(study)

    def signatures(group: object) -> list[str]:
        return [",".join(sorted(space)) for space in group.search_spaces]

    direct_group = optuna.search_space._SearchSpaceGroup()
    direct_group.add_distributions({"x": base["x"], "y": base["fixed"]})
    direct_group.add_distributions(
        {"x": base["x"], "z": distributions.IntDistribution(1, 3)}
    )
    grouped_study = optuna.create_study()
    grouped_study.add_trials(
        [
            frozen(
                TrialState.COMPLETE,
                {"x": base["x"], "fixed": base["fixed"]},
            ),
            frozen(
                TrialState.COMPLETE,
                {"x": base["x"], "z": distributions.IntDistribution(1, 3)},
            ),
        ]
    )
    decomposed = optuna.search_space._GroupDecomposedSearchSpace().calculate(
        grouped_study
    )
    return {
        "without_pruned": list(without_pruned),
        "with_pruned": list(with_pruned),
        "calculator": list(calculated),
        "single_distribution_is_included": "fixed" in without_pruned,
        "group": {
            "calculated_signatures": signatures(decomposed),
            "direct_signatures": signatures(direct_group),
        },
    }


def _enum_contract() -> dict[str, object]:
    def integer_api(item: object) -> dict[str, object]:
        enum_type = type(item)
        value = int(item)
        return {
            "as_integer_ratio": list(item.as_integer_ratio()),
            "bit_count": item.bit_count(),
            "bit_length": item.bit_length(),
            "conjugate": item.conjugate(),
            "denominator": item.denominator,
            "from_bytes_name": enum_type.from_bytes(
                bytes([0, value]), "big"
            ).name,
            "imag": item.imag,
            "is_integer": item.is_integer(),
            "name": item.name,
            "numerator": item.numerator,
            "real": item.real,
            "to_bytes": list(item.to_bytes(2, "big")),
            "value": item.value,
        }

    return {
        "study_direction": [
            {"name": item.name, "value": item.value}
            for item in optuna.study.StudyDirection
        ],
        "trial_state": [
            {
                "name": item.name,
                "value": item.value,
                "is_finished": item.is_finished(),
            }
            for item in TrialState
        ],
        "integer_api": {
            "study_direction": integer_api(optuna.study.StudyDirection.MAXIMIZE),
            "trial_state": integer_api(TrialState.WAITING),
        },
    }


def _unfinished_trial_contract() -> dict[str, object]:
    distribution = {"x": optuna.distributions.FloatDistribution(0.0, 1.0)}
    waiting = optuna.trial.create_trial(
        state=TrialState.WAITING,
        params={"x": 0.25},
        distributions=distribution,
        user_attrs={"source": "waiting"},
    )
    running = optuna.trial.create_trial(
        state=TrialState.RUNNING,
        params={"x": 0.75},
        distributions=distribution,
        user_attrs={"source": "running"},
    )
    study = optuna.create_study(sampler=optuna.samplers.RandomSampler(seed=9))
    study.add_trials([waiting, running])
    before = study.get_trials(deepcopy=False)
    claimed = study.ask()
    claimed_value = claimed.suggest_float("x", 0.0, 1.0)
    after_claim = study.get_trials(deepcopy=False)
    fresh = study.ask()
    return {
        "factory": {
            "waiting_number": waiting.number,
            "waiting_has_start": waiting.datetime_start is not None,
            "waiting_has_complete": waiting.datetime_complete is not None,
            "waiting_has_duration": waiting.duration is not None,
            "running_number": running.number,
            "running_has_start": running.datetime_start is not None,
            "running_has_complete": running.datetime_complete is not None,
            "running_has_duration": running.duration is not None,
        },
        "before": [
            {
                "number": trial.number,
                "state": trial.state.name,
                "params": trial.params,
                "has_start": trial.datetime_start is not None,
                "has_complete": trial.datetime_complete is not None,
                "has_duration": trial.duration is not None,
            }
            for trial in before
        ],
        "claimed_number": claimed.number,
        "claimed_value": float(claimed_value),
        "claimed_params": claimed.params,
        "claimed_user_attrs": claimed.user_attrs,
        "after_claim_states": [trial.state.name for trial in after_claim],
        "after_claim_has_start": [
            trial.datetime_start is not None for trial in after_claim
        ],
        "fresh_number": fresh.number,
        "fresh_params": fresh.params,
    }


def _base_trial_contract() -> dict[str, object]:
    base = optuna.trial.BaseTrial
    fixed = optuna.trial.FixedTrial({"x": 0.5})
    frozen = optuna.trial.create_trial(value=1.0)
    live = optuna.create_study().ask()
    try:
        base()
    except Exception as error:  # noqa: BLE001
        construction_error = type(error).__name__
    else:
        construction_error = None
    return {
        "construction_error": construction_error,
        "is_base_trial": {
            "fixed": isinstance(fixed, base),
            "frozen": isinstance(frozen, base),
            "trial": isinstance(live, base),
        },
        "numbers": {
            "fixed": fixed.number,
            "frozen": frozen.number,
            "trial": live.number,
        },
    }


def _base_component_contract() -> dict[str, object]:
    sampler = optuna.samplers.RandomSampler(seed=5)
    pruner = optuna.pruners.NopPruner()
    study = optuna.create_study()
    frozen = optuna.trial.create_trial(value=1.0)
    errors: dict[str, str | None] = {}
    for name, operation in {
        "sampler": optuna.samplers.BaseSampler,
        "pruner": optuna.pruners.BasePruner,
    }.items():
        try:
            operation()
        except Exception as error:  # noqa: BLE001
            errors[name] = type(error).__name__
        else:
            errors[name] = None
    return {
        "sampler_is_base": isinstance(sampler, optuna.samplers.BaseSampler),
        "pruner_is_base": isinstance(pruner, optuna.pruners.BasePruner),
        "nop_decision": pruner.prune(study, frozen),
        "construction_errors": errors,
    }


def _qmc_warning_contract() -> dict[str, object]:
    logger = logging.getLogger("optuna.samplers._qmc")

    class Records(logging.Handler):
        def __init__(self) -> None:
            super().__init__()
            self.messages: list[str] = []

        def emit(self, record: logging.LogRecord) -> None:
            self.messages.append(record.getMessage())

    def capture(operation: object) -> list[str]:
        handler = Records()
        logger.addHandler(handler)
        try:
            operation()  # type: ignore[operator]
        finally:
            logger.removeHandler(handler)
        return handler.messages

    asynchronous_enabled = capture(
        lambda: optuna.samplers.QMCSampler(
            qmc_type="sobol",
            scramble=True,
            seed=None,
            warn_asynchronous_seeding=True,
        )
    )
    asynchronous_disabled = capture(
        lambda: optuna.samplers.QMCSampler(
            qmc_type="sobol",
            scramble=True,
            seed=None,
            warn_asynchronous_seeding=False,
        )
    )

    def sample_categorical(warn: bool) -> None:
        sampler = optuna.samplers.QMCSampler(
            seed=11, warn_independent_sampling=warn
        )
        study = optuna.create_study(sampler=sampler)
        for _ in range(2):
            trial = study.ask()
            value = trial.suggest_categorical("kind", ["a", "b"])
            study.tell(trial, 0.0 if value == "a" else 1.0)

    independent_enabled = capture(lambda: sample_categorical(True))
    independent_disabled = capture(lambda: sample_categorical(False))
    return {
        "asynchronous_enabled_count": len(asynchronous_enabled),
        "asynchronous_disabled_count": len(asynchronous_disabled),
        "independent_enabled_count": len(independent_enabled),
        "independent_disabled_count": len(independent_disabled),
    }


def _distribution_json_contract() -> dict[str, object]:
    distributions = optuna.distributions
    current = {
        "float": distributions.FloatDistribution(0.0, 1.0),
        "log_float": distributions.FloatDistribution(0.001, 1.0, log=True),
        "stepped_float": distributions.FloatDistribution(0.0, 1.0, step=0.2),
        "integer": distributions.IntDistribution(1, 9, step=2),
        "log_integer": distributions.IntDistribution(1, 100, log=True),
        "categorical": distributions.CategoricalDistribution(["A", 2, True]),
    }
    aliases: dict[str, object] = {}
    warning_categories: dict[str, list[str]] = {}
    alias_specs = {
        "uniform": (distributions.UniformDistribution, (0.0, 1.0)),
        "log_uniform": (distributions.LogUniformDistribution, (0.001, 1.0)),
        "discrete_uniform": (
            distributions.DiscreteUniformDistribution,
            (0.0, 1.0, 0.2),
        ),
        "int_uniform": (distributions.IntUniformDistribution, (1, 9, 2)),
        "int_log_uniform": (distributions.IntLogUniformDistribution, (1, 100, 1)),
    }
    for name, (constructor, args) in alias_specs.items():
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            aliases[name] = constructor(*args)
        warning_categories[name] = [type(item.message).__name__ for item in caught]

    errors: dict[str, str] = {}
    pairs = {
        "kind": (
            distributions.FloatDistribution(0.0, 1.0),
            distributions.IntDistribution(0, 1),
        ),
        "log": (
            distributions.FloatDistribution(0.0, 1.0),
            distributions.FloatDistribution(0.001, 1.0, log=True),
        ),
        "categorical": (
            distributions.CategoricalDistribution(["A", "B"]),
            distributions.CategoricalDistribution(["A", "C"]),
        ),
    }
    for name, pair in pairs.items():
        try:
            distributions.check_distribution_compatibility(*pair)
        except Exception as error:  # noqa: BLE001 - record public exception type.
            errors[name] = type(error).__name__
    distributions.check_distribution_compatibility(
        distributions.FloatDistribution(0.0, 1.0),
        distributions.FloatDistribution(-1.0, 2.0, step=0.2),
    )
    all_distributions = {**current, **aliases}
    encoded = {
        name: distributions.distribution_to_json(value)
        for name, value in all_distributions.items()
    }
    return {
        "encoded": encoded,
        "roundtrip_types": {
            name: type(distributions.json_to_distribution(value)).__name__
            for name, value in encoded.items()
        },
        "alias_warning_categories": warning_categories,
        "compatibility": {
            "range_change_allowed": True,
            "errors": errors,
        },
    }


def _grid_trials() -> list[dict[str, object]]:
    sampler = optuna.samplers.GridSampler(
        {"x": [-1.0, 0.0, 1.0], "mode": ["A", "B"]}, seed=17
    )
    study = optuna.create_study(sampler=sampler)
    rows: list[dict[str, object]] = []

    def objective(trial: optuna.Trial) -> float:
        row = {
            "x": trial.suggest_float("x", -1.0, 1.0),
            "mode": trial.suggest_categorical("mode", ["A", "B"]),
        }
        rows.append(row)
        return float(row["x"])

    study.optimize(objective, n_trials=6)
    return rows


def _nsgaii_trials() -> list[dict[str, object]]:
    sampler = optuna.samplers.NSGAIISampler(seed=19, population_size=4)
    study = optuna.create_study(
        directions=["minimize", "minimize"], sampler=sampler
    )
    rows: list[dict[str, object]] = []

    def objective(trial: optuna.Trial) -> list[float]:
        x = trial.suggest_float("x", -1.0, 1.0)
        mesh = trial.suggest_int("mesh", 1, 5, step=2)
        mode = trial.suggest_categorical("mode", ["A", "B"])
        rows.append({"x": x, "mesh": mesh, "mode": mode})
        return [x * x + 0.1 * mesh, (x - 0.5) ** 2 + (0.2 if mode == "B" else 0.0)]

    study.optimize(objective, n_trials=32)
    return rows


def _nsgaiii_trials() -> list[dict[str, object]]:
    sampler = optuna.samplers.NSGAIIISampler(seed=23, population_size=4)
    study = optuna.create_study(
        directions=["minimize", "minimize", "minimize"], sampler=sampler
    )
    rows: list[dict[str, object]] = []

    def objective(trial: optuna.Trial) -> list[float]:
        x = trial.suggest_float("x", -1.0, 1.0)
        mesh = trial.suggest_int("mesh", 1, 5, step=2)
        mode = trial.suggest_categorical("mode", ["A", "B"])
        rows.append({"x": x, "mesh": mesh, "mode": mode})
        return [
            x * x + 0.1 * mesh,
            (x - 0.5) ** 2 + (0.2 if mode == "B" else 0.0),
            (x + 0.25) ** 2 + 0.05 * mesh,
        ]

    study.optimize(objective, n_trials=32)
    return rows


def _brute_force_trials() -> list[dict[str, object]]:
    sampler = optuna.samplers.BruteForceSampler(seed=29)
    study = optuna.create_study(sampler=sampler)
    rows: list[dict[str, object]] = []

    def objective(trial: optuna.Trial) -> float:
        mesh = trial.suggest_int("mesh", 1, 3)
        mode = trial.suggest_categorical("mode", ["A", "B"])
        rows.append({"mesh": mesh, "mode": mode})
        return float(mesh) + (0.1 if mode == "B" else 0.0)

    study.optimize(objective)
    return rows


def _conditional_brute_force_trials() -> list[dict[str, str]]:
    sampler = optuna.samplers.BruteForceSampler(seed=79)
    study = optuna.create_study(sampler=sampler)
    rows: list[dict[str, str]] = []

    def objective(trial: optuna.Trial) -> float:
        branch = trial.suggest_categorical("branch", ["depth", "mode"])
        if branch == "depth":
            value = trial.suggest_int("depth", 1, 2)
            rows.append({"branch": branch, "parameter": "depth", "value": str(value)})
            return float(value)
        value = trial.suggest_categorical("mode", ["A", "B", "C"])
        rows.append({"branch": branch, "parameter": "mode", "value": value})
        return float(ord(value) - ord("A"))

    study.optimize(objective)
    return rows


def _cmaes_trials() -> list[dict[str, float]]:
    sampler = optuna.samplers.CmaEsSampler(
        seed=31, n_startup_trials=1, popsize=4
    )
    study = optuna.create_study(sampler=sampler)
    rows: list[dict[str, float]] = []

    def objective(trial: optuna.Trial) -> float:
        x = trial.suggest_float("x", -2.0, 2.0)
        y = trial.suggest_float("y", -1.0, 3.0)
        rows.append({"x": x, "y": y})
        return (x - 0.4) ** 2 + 0.5 * (y + 0.2) ** 2

    study.optimize(objective, n_trials=32)
    return rows


def _cmaes_independent_sampler_trials() -> list[dict[str, object]]:
    sampler = optuna.samplers.CmaEsSampler(
        seed=31,
        n_startup_trials=1,
        popsize=4,
        independent_sampler=optuna.samplers.RandomSampler(seed=211),
        warn_independent_sampling=False,
    )
    study = optuna.create_study(sampler=sampler)
    rows: list[dict[str, object]] = []
    for _ in range(24):
        trial = study.ask()
        x = trial.suggest_float("x", -2.0, 2.0)
        y = trial.suggest_float("y", -1.0, 3.0)
        mode = trial.suggest_categorical("mode", ["A", "B", "C"])
        loss = (x - 0.4) ** 2 + 0.5 * (y + 0.2) ** 2 + 0.1 * "ABC".index(mode)
        study.tell(trial, loss)
        rows.append({"x": x, "y": y, "mode": mode})
    return rows


def _scrambled_qmc_trials() -> dict[str, list[dict[str, float]]]:
    result: dict[str, list[dict[str, float]]] = {}
    for qmc_type in ("sobol", "halton"):
        sampler = optuna.samplers.QMCSampler(
            qmc_type=qmc_type, scramble=True, seed=47
        )
        study = optuna.create_study(sampler=sampler)
        rows: list[dict[str, float]] = []
        for _ in range(16):
            trial = study.ask()
            x = trial.suggest_float("x", -1.0, 1.0)
            y = trial.suggest_float("y", 0.0, 4.0)
            study.tell(trial, x * x + y * y)
            rows.append({"x": x, "y": y})
        result[qmc_type] = rows
    return result


def _gp_trials() -> list[dict[str, object]]:
    sampler = optuna.samplers.GPSampler(seed=53, n_startup_trials=10)
    study = optuna.create_study(sampler=sampler)
    rows: list[dict[str, object]] = []
    for _ in range(12):
        trial = study.ask()
        x = trial.suggest_float("x", -1.0, 1.0)
        mesh = trial.suggest_int("mesh", 1, 5, step=2)
        mode = trial.suggest_categorical("mode", ["A", "B"])
        study.tell(trial, x * x + 0.1 * mesh + (0.2 if mode == "B" else 0.0))
        rows.append({"x": x, "mesh": mesh, "mode": mode})
    return rows


def _gp_constraint_trials() -> list[dict[str, float]]:
    sampler = optuna.samplers.GPSampler(
        seed=89,
        n_startup_trials=5,
        constraints_func=lambda trial: trial.user_attrs["constraints"],
    )
    study = optuna.create_study(sampler=sampler)
    rows: list[dict[str, float]] = []
    for _ in range(8):
        trial = study.ask()
        x = trial.suggest_float("x", -1.0, 1.0)
        constraint = x - 0.1
        trial.set_user_attr("constraints", [constraint])
        study.tell(trial, (x - 0.35) ** 2)
        rows.append({"x": x, "constraint": constraint})
    return rows


def _nsgaii_crossover_trials() -> dict[str, list[dict[str, float]]]:
    crossovers = {
        "uniform": optuna.samplers.nsgaii.UniformCrossover(),
        "blxalpha": optuna.samplers.nsgaii.BLXAlphaCrossover(),
        "sbx": optuna.samplers.nsgaii.SBXCrossover(),
        "vsbx": optuna.samplers.nsgaii.VSBXCrossover(),
        "spx": optuna.samplers.nsgaii.SPXCrossover(),
        "undx": optuna.samplers.nsgaii.UNDXCrossover(),
    }
    result: dict[str, list[dict[str, float]]] = {}
    for name, crossover in crossovers.items():
        sampler = optuna.samplers.NSGAIISampler(
            seed=73,
            population_size=4,
            mutation_prob=0.0,
            crossover_prob=1.0,
            crossover=crossover,
        )
        study = optuna.create_study(
            sampler=sampler, directions=["minimize", "minimize"]
        )
        rows: list[dict[str, float]] = []
        for _ in range(12):
            trial = study.ask()
            x = trial.suggest_float("x", -1.0, 1.0)
            y = trial.suggest_float("y", -2.0, 2.0)
            z = trial.suggest_float("z", 0.0, 3.0)
            study.tell(
                trial,
                [x * x + 0.2 * y * y + 0.1 * z, (x - 0.4) ** 2 + (y + 0.3) ** 2 + z * z],
            )
            rows.append({"x": x, "y": y, "z": z})
        result[name] = rows
    return result


def _partial_fixed_trials() -> list[dict[str, float]]:
    sampler = optuna.samplers.PartialFixedSampler(
        {"x": 0.25}, optuna.samplers.RandomSampler(seed=61)
    )
    study = optuna.create_study(sampler=sampler)
    rows: list[dict[str, float]] = []
    for _ in range(5):
        trial = study.ask()
        x = trial.suggest_float("x", 0.0, 1.0)
        y = trial.suggest_float("y", -1.0, 1.0)
        study.tell(trial, x * x + y * y)
        rows.append({"x": x, "y": y})
    return rows


def _multivariate_tpe_trials() -> list[dict[str, object]]:
    sampler = optuna.samplers.TPESampler(
        seed=67, n_startup_trials=4, multivariate=True
    )
    study = optuna.create_study(sampler=sampler)
    rows: list[dict[str, object]] = []
    for _ in range(30):
        trial = study.ask()
        x = trial.suggest_float("x", -1.0, 1.0)
        mesh = trial.suggest_int("mesh", 1, 5, step=2)
        mode = trial.suggest_categorical("mode", ["A", "B"])
        study.tell(
            trial,
            (x - 0.2) ** 2
            + 0.05 * mesh
            + (0.0 if mode == "B" else 0.2),
        )
        rows.append({"x": x, "mesh": mesh, "mode": mode})
    return rows


def _unscrambled_qmc_trials() -> dict[str, list[dict[str, float]]]:
    result: dict[str, list[dict[str, float]]] = {}
    for qmc_type in ("sobol", "halton"):
        sampler = optuna.samplers.QMCSampler(
            qmc_type=qmc_type, scramble=False, seed=71
        )
        study = optuna.create_study(sampler=sampler)
        rows: list[dict[str, float]] = []
        for _ in range(16):
            trial = study.ask()
            x = trial.suggest_float("x", -1.0, 1.0)
            y = trial.suggest_float("y", 0.0, 4.0)
            study.tell(trial, x * x + y * y)
            rows.append({"x": x, "y": y})
        result[qmc_type] = rows
    return result


def _pruner_contract() -> dict[str, object]:
    def add_completed(
        study: optuna.Study, steps: list[int], values: list[float]
    ) -> None:
        trial = study.ask()
        for step, value in zip(steps, values, strict=True):
            trial.report(value, step)
        study.tell(trial, values[-1])

    percentile = optuna.create_study(
        pruner=optuna.pruners.PercentilePruner(
            50.0,
            n_startup_trials=0,
            n_warmup_steps=0,
            interval_steps=2,
            n_min_trials=1,
        )
    )
    add_completed(percentile, [1, 3], [1.0, 1.0])
    add_completed(percentile, [1, 3], [3.0, 3.0])
    percentile_trial = percentile.ask()
    percentile_trial.report(5.0, 1)
    percentile_trial.report(4.0, 3)

    maximize = optuna.create_study(
        direction="maximize",
        pruner=optuna.pruners.PercentilePruner(50.0, n_startup_trials=0),
    )
    add_completed(maximize, [0], [1.0])
    add_completed(maximize, [0], [3.0])
    maximize_trial = maximize.ask()
    maximize_trial.report(1.5, 0)

    median = optuna.create_study(
        pruner=optuna.pruners.MedianPruner(n_startup_trials=0)
    )
    add_completed(median, [0], [1.0])
    add_completed(median, [0], [3.0])
    median_trial = median.ask()
    median_trial.report(5.0, 0)

    threshold = optuna.create_study(
        pruner=optuna.pruners.ThresholdPruner(lower=0.0, upper=10.0)
    )
    threshold_trial = threshold.ask()
    threshold_trial.report(5.0, 0)
    threshold_first = threshold_trial.should_prune()
    threshold_trial.report(11.0, 1)
    threshold_second = threshold_trial.should_prune()
    threshold_nan = threshold.ask()
    threshold_nan.report(float("nan"), 0)

    patient = optuna.create_study(
        pruner=optuna.pruners.PatientPruner(None, patience=1, min_delta=0.0)
    )
    patient_trial = patient.ask()
    for step, value in enumerate([10.0, 4.0, 5.0, 6.0]):
        patient_trial.report(value, step)
    wrapped = optuna.create_study(
        pruner=optuna.pruners.PatientPruner(
            optuna.pruners.NopPruner(), patience=1
        )
    )
    wrapped_trial = wrapped.ask()
    for step, value in enumerate([10.0, 4.0, 5.0, 6.0]):
        wrapped_trial.report(value, step)

    halving = optuna.create_study(
        pruner=optuna.pruners.SuccessiveHalvingPruner(
            min_resource=1, reduction_factor=2
        )
    )
    first = halving.ask()
    first.report(1.0, 1)
    first_decision = first.should_prune()
    first_rung = first.system_attrs["completed_rung_0"]
    halving.tell(first, 1.0)
    second = halving.ask()
    second.report(2.0, 1)
    second_decision = second.should_prune()
    second_rung = second.system_attrs["completed_rung_0"]

    bootstrap = optuna.create_study(
        pruner=optuna.pruners.SuccessiveHalvingPruner(
            min_resource=1, reduction_factor=2, bootstrap_count=1
        )
    )
    bootstrap_trial = bootstrap.ask()
    bootstrap_trial.report(1.0, 1)

    hyperband_pruner = optuna.pruners.HyperbandPruner(
        min_resource=1, max_resource=9, reduction_factor=3
    )
    hyperband = optuna.create_study(study_name="hb", pruner=hyperband_pruner)
    hyperband_trial = hyperband.ask()
    hyperband_trial.report(1.0, 0)
    hyperband_decision = hyperband_trial.should_prune()
    bracket_ids = [
        hyperband_pruner._get_bracket_id(hyperband, trial)
        for trial in [
            optuna.trial.FrozenTrial(
                number=number,
                state=TrialState.RUNNING,
                value=None,
                datetime_start=None,
                datetime_complete=None,
                params={},
                distributions={},
                user_attrs={},
                system_attrs={},
                intermediate_values={},
                trial_id=number,
                values=None,
            )
            for number in range(10)
        ]
    ]

    wilcoxon = optuna.create_study(
        pruner=optuna.pruners.WilcoxonPruner(p_threshold=0.1, n_startup_steps=2)
    )
    add_completed(wilcoxon, list(range(6)), [0.0] * 6)
    wilcoxon_trial = wilcoxon.ask()
    for step in range(6):
        wilcoxon_trial.report(10.0, step)
    nonfinite = wilcoxon.ask()
    nonfinite.report(float("inf"), 0)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        nonfinite_decision = nonfinite.should_prune()

    return {
        "percentile_minimize": bool(percentile_trial.should_prune()),
        "percentile_maximize": bool(maximize_trial.should_prune()),
        "median": bool(median_trial.should_prune()),
        "threshold": [
            bool(threshold_first),
            bool(threshold_second),
            bool(threshold_nan.should_prune()),
        ],
        "patient": [
            bool(patient_trial.should_prune()),
            bool(wrapped_trial.should_prune()),
        ],
        "successive_halving": {
            "decisions": [bool(first_decision), bool(second_decision)],
            "rung_values": [float(first_rung), float(second_rung)],
            "bootstrap": bool(bootstrap_trial.should_prune()),
        },
        "hyperband": {
            "first_decision": bool(hyperband_decision),
            "bracket_ids": [int(value) for value in bracket_ids],
        },
        "wilcoxon": [
            bool(wilcoxon_trial.should_prune()),
            bool(nonfinite_decision),
        ],
    }


def _constraint_contract() -> dict[str, object]:
    sampler = optuna.samplers.NSGAIISampler(
        seed=73, population_size=4, constraints_func=lambda trial: trial.user_attrs["c"]
    )
    study = optuna.create_study(
        directions=["minimize", "minimize"], sampler=sampler
    )
    values = [[0.0, 0.0], [1.0, 2.0], [2.0, 1.0], [-1.0, -1.0]]
    constraints = [[1.0], [-1.0], [0.0], [2.0]]
    for objective_values, constraint_values in zip(values, constraints, strict=True):
        trial = study.ask()
        trial.set_user_attr("c", constraint_values)
        study.tell(trial, objective_values)
    return {
        "pareto_trial_numbers": sorted(trial.number for trial in study.best_trials),
        "states": [trial.state.name for trial in study.trials],
        "constraints": [list(trial.system_attrs["constraints"]) for trial in study.trials],
    }


def _tell_contract() -> dict[str, object]:
    study = optuna.create_study()
    pruned = study.ask()
    pruned.report(3.5, 0)
    pruned.report(2.5, 1)
    frozen_pruned = study.tell(pruned, state=TrialState.PRUNED)
    failed = study.ask()
    frozen_failed = study.tell(failed, state=TrialState.FAIL)
    missing = study.ask()
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        frozen_missing = study.tell(missing)
    infinite = study.ask()
    frozen_infinite = study.tell(infinite, float("inf"))
    by_number = study.ask()
    by_number.report(4.5, 1)
    frozen_by_number = study.tell(by_number.number, 4.25)
    return {
        "pruned_state": frozen_pruned.state.name,
        "pruned_value": frozen_pruned.value,
        "failed_state": frozen_failed.state.name,
        "failed_value_is_none": frozen_failed.value is None,
        "missing_state": frozen_missing.state.name,
        "infinite_state": frozen_infinite.state.name,
        "infinite_value": "Infinity" if frozen_infinite.value == float("inf") else None,
        "by_number": {
            "number": frozen_by_number.number,
            "state": frozen_by_number.state.name,
            "value": frozen_by_number.value,
            "last_step": frozen_by_number.last_step,
        },
    }


def _trial_pruned_exception_contract() -> dict[str, object]:
    callback_rows: list[dict[str, object]] = []

    def objective(trial: optuna.Trial) -> float:
        trial.report(7.0, 0)
        trial.report(3.0, 2)
        raise optuna.TrialPruned("oracle prune")

    def callback(_: optuna.Study, trial: optuna.trial.FrozenTrial) -> None:
        callback_rows.append(
            {
                "state": trial.state.name,
                "value": trial.value,
                "last_step": trial.last_step,
            }
        )

    study = optuna.create_study()
    study.optimize(objective, n_trials=1, callbacks=[callback])
    trial = study.trials[0]
    return {
        "state": trial.state.name,
        "value": trial.value,
        "last_step": trial.last_step,
        "callback_count": len(callback_rows),
        "callback": callback_rows[0],
    }


def _lifecycle_error_contract() -> dict[str, object]:
    finished = optuna.create_study()
    trial = finished.ask()
    finished.tell(trial, 1.0)
    repeat_error = ""
    try:
        finished.tell(trial, 2.0)
    except Exception as error:  # noqa: BLE001
        repeat_error = type(error).__name__
    skipped = finished.tell(trial, 999.0, skip_if_finished=True)

    def failing_constraint(_: optuna.trial.FrozenTrial) -> list[float]:
        raise RuntimeError("constraint callback failed")

    constrained = optuna.create_study(
        sampler=optuna.samplers.TPESampler(seed=83, constraints_func=failing_constraint)
    )
    constraint_error = ""
    try:
        constrained.optimize(lambda item: item.suggest_float("x", 0.0, 1.0), n_trials=1)
    except Exception as error:  # noqa: BLE001
        constraint_error = type(error).__name__
    constrained_trial = constrained.trials[0]
    return {
        "finished_tell": {
            "repeat_error": repeat_error,
            "skip_value": float(skipped.value),
            "skip_state": skipped.state.name,
        },
        "constraint_callback_failure": {
            "error": constraint_error,
            "state": constrained_trial.state.name,
            "value_is_finite": bool(np.isfinite(constrained_trial.value)),
            "constraint_is_none": constrained_trial.system_attrs["constraints"] is None,
        },
    }


def _fixed_trial_contract() -> dict[str, object]:
    trial = optuna.trial.FixedTrial(
        {"x": 0.5, "n": 3, "kind": "b"}, number=7
    )
    initial_params = dict(trial.params)
    values = {
        "x": float(trial.suggest_float("x", 0.0, 1.0)),
        "n": int(trial.suggest_int("n", 1, 5)),
        "kind": str(trial.suggest_categorical("kind", ["a", "b"])),
    }
    trial.report(2.5, 3)
    trial.set_user_attr("owner", "matlab")
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", FutureWarning)
        trial.set_system_attr("generation", 2)
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        repeated = float(trial.suggest_float("x", 0.6, 1.0))

    errors: dict[str, str] = {}
    for name, operation in {
        "missing": lambda: trial.suggest_float("missing", 0.0, 1.0),
        "different_kind": lambda: trial.suggest_int("x", 0, 2),
        "categorical_choice": lambda: optuna.trial.FixedTrial(
            {"choice": "z"}
        ).suggest_categorical("choice", ["a", "b"]),
    }.items():
        try:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                operation()
        except Exception as error:  # noqa: BLE001
            errors[name] = type(error).__name__
    return {
        "number": trial.number,
        "initial_params": initial_params,
        "initial_distribution_count": 0,
        "values": values,
        "should_prune": bool(trial.should_prune()),
        "repeated_out_of_range": repeated,
        "warning_categories": [item.category.__name__ for item in caught],
        "distribution_types": {
            name: type(distribution).__name__
            for name, distribution in trial.distributions.items()
        },
        "user_attrs": dict(trial.user_attrs),
        "system_attrs": dict(trial.system_attrs),
        "errors": errors,
    }


def _frozen_trial_contract() -> dict[str, object]:
    trial = optuna.trial.create_trial(
        value=1.2,
        params={"x": 0.5, "n": 3, "kind": "b"},
        distributions={
            "x": optuna.distributions.FloatDistribution(0.0, 1.0),
            "n": optuna.distributions.IntDistribution(1, 5),
            "kind": optuna.distributions.CategoricalDistribution(["a", "b"]),
        },
        user_attrs={"owner": "upstream"},
        system_attrs={"generation": 1},
        intermediate_values={2: 3.0},
    )
    values = {
        "x": float(trial.suggest_float("x", 0.0, 1.0)),
        "n": int(trial.suggest_int("n", 1, 5)),
        "kind": str(trial.suggest_categorical("kind", ["a", "b"])),
    }
    trial.report(99.0, 9)
    trial.set_user_attr("owner", "matlab")
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", FutureWarning)
        trial.set_system_attr("generation", 2)
    errors: dict[str, str] = {}
    for name, operation in {
        "missing": lambda: trial.suggest_float("missing", 0.0, 1.0),
        "different_kind": lambda: trial.suggest_int("x", 0, 2),
        "categorical_choices": lambda: trial.suggest_categorical(
            "kind", ["a", "c"]
        ),
    }.items():
        try:
            operation()
        except Exception as error:  # noqa: BLE001
            errors[name] = type(error).__name__
    return {
        "values": values,
        "last_step": trial.last_step,
        "should_prune": trial.should_prune(),
        "report_is_noop": 9 not in trial.intermediate_values,
        "user_owner": trial.user_attrs["owner"],
        "system_generation": trial.system_attrs["generation"],
        "errors": errors,
    }


def _study_management_contract() -> dict[str, object]:
    source_storage = optuna.storages.InMemoryStorage()
    target_storage = optuna.storages.InMemoryStorage()
    source = optuna.create_study(
        study_name="alpha", direction="maximize", storage=source_storage
    )
    source.add_trial(
        optuna.trial.create_trial(
            value=1.25,
            params={"x": 0.5},
            distributions={
                "x": optuna.distributions.FloatDistribution(0.0, 1.0)
            },
            user_attrs={"origin": "oracle"},
        )
    )
    loaded = optuna.load_study(study_name="alpha", storage=source_storage)
    loaded_result = {
        "name": loaded.study_name,
        "direction": loaded.direction.name,
        "trial_count": len(loaded.trials),
        "best_value": float(loaded.best_value),
    }
    optuna.copy_study(
        from_study_name="alpha",
        from_storage=source_storage,
        to_storage=target_storage,
        to_study_name="beta",
    )
    copied = optuna.load_study(study_name="beta", storage=target_storage)
    summary = optuna.get_all_study_summaries(target_storage)[0]
    summary_without_best = optuna.get_all_study_summaries(
        target_storage, include_best_trial=False
    )[0]
    multi_storage = optuna.storages.InMemoryStorage()
    optuna.create_study(
        study_name="multi",
        directions=["minimize", "maximize"],
        storage=multi_storage,
    )
    multi_summary = optuna.get_all_study_summaries(multi_storage)[0]
    try:
        _ = multi_summary.direction
    except Exception as error:  # noqa: BLE001
        multi_direction_error = type(error).__name__
    else:
        multi_direction_error = None
    names_before_delete = optuna.get_all_study_names(source_storage)
    optuna.delete_study(study_name="alpha", storage=source_storage)
    names_after_delete = optuna.get_all_study_names(source_storage)
    return {
        "loaded": loaded_result,
        "copied": {
            "name": copied.study_name,
            "direction": copied.direction.name,
            "trial_count": len(copied.trials),
            "best_value": float(copied.best_value),
            "param_x": float(copied.best_trial.params["x"]),
            "user_origin": str(copied.best_trial.user_attrs["origin"]),
        },
        "summary": {
            "type": type(summary).__name__,
            "name": summary.study_name,
            "direction": summary.direction.name,
            "directions": [item.name for item in summary.directions],
            "trial_count": summary.n_trials,
            "best_value": float(summary.best_trial.value),
            "without_best_has_best": summary_without_best.best_trial is not None,
        },
        "multi_summary": {
            "directions": [item.name for item in multi_summary.directions],
            "direction_error": multi_direction_error,
        },
        "names_before_delete": names_before_delete,
        "names_after_delete": names_after_delete,
    }


def _distribution_public_member_contract() -> dict[str, object]:
    distributions = optuna.distributions
    try:
        distributions.BaseDistribution()
    except Exception as error:  # noqa: BLE001 - abstract public API contract.
        base_error = type(error).__name__
    else:
        raise RuntimeError("Optuna BaseDistribution unexpectedly became concrete.")

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        instances = {
            "FloatDistribution": distributions.FloatDistribution(1.0, 1.0),
            "IntDistribution": distributions.IntDistribution(1, 1),
            "CategoricalDistribution": distributions.CategoricalDistribution(
                ["A", "B", 3]
            ),
            "UniformDistribution": distributions.UniformDistribution(1.0, 1.0),
            "LogUniformDistribution": distributions.LogUniformDistribution(1.0, 2.0),
            "DiscreteUniformDistribution": distributions.DiscreteUniformDistribution(
                0.0, 1.0, 0.2
            ),
            "IntUniformDistribution": distributions.IntUniformDistribution(1, 3, 1),
            "IntLogUniformDistribution": distributions.IntLogUniformDistribution(
                1, 3, 1
            ),
        }

    single = {name: value.single() for name, value in instances.items()}
    registry = [value.__name__ for value in distributions.DISTRIBUTION_CLASSES]
    return {
        "base_construction_error": base_error,
        "categorical_choice_type": str(distributions.CategoricalChoiceType),
        "categorical_external": instances[
            "CategoricalDistribution"
        ].to_external_repr(2),
        "categorical_internal": instances[
            "CategoricalDistribution"
        ].to_internal_repr("B"),
        "discrete_q": instances["DiscreteUniformDistribution"].q,
        "distribution_classes": registry,
        "float_external": instances["FloatDistribution"].to_external_repr(1.25),
        "float_internal": instances["FloatDistribution"].to_internal_repr("1.25"),
        "int_external": instances["IntDistribution"].to_external_repr(3.9),
        "int_internal": instances["IntDistribution"].to_internal_repr("3"),
        "single": single,
    }


def _normalize_dataframe_value(value: object) -> object:
    if isinstance(value, TrialState):
        return value.name
    if isinstance(value, np.generic):
        value = value.item()
    if isinstance(value, float) and np.isnan(value):
        return None
    return value


def _dataframe_snapshot(frame: object) -> dict[str, object]:
    columns = list(frame.columns)  # type: ignore[attr-defined]
    flat_columns = [
        "_".join(str(part) for part in column if str(part))
        if isinstance(column, tuple)
        else str(column)
        for column in columns
    ]
    values = {
        flat_columns[index]: [
            _normalize_dataframe_value(value)
            for value in frame.iloc[:, index].tolist()  # type: ignore[attr-defined]
        ]
        for index in range(len(columns))
    }
    return {
        "flat_columns": flat_columns,
        "multi_columns": [
            {
                "top": str(column[0]),
                "sub": str(column[1]),
            }
            if isinstance(column, tuple)
            else {"top": str(column), "sub": ""}
            for column in columns
        ],
        "values": values,
    }


def _trials_dataframe_contract() -> dict[str, object]:
    distributions = optuna.distributions
    common_distributions = {
        "x": distributions.FloatDistribution(0.0, 1.0),
        "mode": distributions.CategoricalDistribution(["A", "B"]),
    }
    trials = [
        optuna.trial.create_trial(
            value=1.25,
            params={"x": 0.25, "mode": "A"},
            distributions=common_distributions,
            user_attrs={"owner": "lab"},
            system_attrs={"origin": "oracle"},
            state=TrialState.COMPLETE,
        ),
        optuna.trial.create_trial(
            params={"x": 0.75, "mode": "B"},
            distributions=common_distributions,
            intermediate_values={0: 3.0, 2: 1.0},
            user_attrs={"owner": "mdx"},
            system_attrs={"origin": "oracle"},
            state=TrialState.PRUNED,
        ),
        optuna.trial.create_trial(
            params={"x": 0.5, "mode": "A"},
            distributions=common_distributions,
            user_attrs={"owner": "lab"},
            system_attrs={"origin": "imported"},
            state=TrialState.FAIL,
        ),
    ]
    study = optuna.create_study()
    study.add_trials(trials)
    attrs = (
        "number",
        "value",
        "params",
        "user_attrs",
        "system_attrs",
        "state",
    )
    flat = study.trials_dataframe(attrs=attrs, multi_index=False)
    multi = study.trials_dataframe(attrs=attrs, multi_index=True)
    single = _dataframe_snapshot(flat)
    single["multi_columns"] = _dataframe_snapshot(multi)["multi_columns"]

    metric_study = optuna.create_study(directions=["minimize", "maximize"])
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        metric_study.set_metric_names(["loss", "gain"])
    metric_study.add_trial(
        optuna.trial.create_trial(
            values=[1.5, 2.5],
            params={"x": 0.5},
            distributions={"x": distributions.FloatDistribution(0.0, 1.0)},
            state=TrialState.COMPLETE,
        )
    )
    metric_attrs = ("number", "value", "params", "state")
    metric_flat = metric_study.trials_dataframe(
        attrs=metric_attrs, multi_index=False
    )
    metric_multi = metric_study.trials_dataframe(
        attrs=metric_attrs, multi_index=True
    )
    metric = _dataframe_snapshot(metric_flat)
    metric["multi_columns"] = _dataframe_snapshot(metric_multi)["multi_columns"]

    single_metric_study = optuna.create_study()
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        single_metric_study.set_metric_names(["loss"])
    single_metric_study.add_trial(optuna.trial.create_trial(value=1.2))
    single_metric_flat = single_metric_study.trials_dataframe(
        attrs=("number", "value", "state"), multi_index=False
    )
    single_metric_multi = single_metric_study.trials_dataframe(
        attrs=("number", "value", "state"), multi_index=True
    )
    single_metric = _dataframe_snapshot(single_metric_flat)
    single_metric["multi_columns"] = _dataframe_snapshot(single_metric_multi)[
        "multi_columns"
    ]

    errors: dict[str, str] = {}
    for name, candidate_attrs in {
        "unknown_attr": ("not_a_frozen_trial_field",),
        "empty_attrs": (),
    }.items():
        try:
            study.trials_dataframe(attrs=candidate_attrs)
        except Exception as error:  # noqa: BLE001 - public error-type oracle.
            errors[name] = type(error).__name__

    empty = optuna.create_study().trials_dataframe(attrs=attrs)
    return {
        "attrs": list(attrs),
        "multi_index_default": False,
        "single": single,
        "single_metric_name": single_metric,
        "metric_names": metric,
        "empty_column_count": len(empty.columns),
        "errors": errors,
    }


def _importance_contract() -> dict[str, object]:
    study = optuna.create_study()
    rows: list[dict[str, object]] = []
    modes = ["A", "B", "C"]
    for index in range(24):
        x = -1.0 + 2.0 * index / 23.0
        y = float((index * 7) % 11) / 10.0
        mode = modes[index % len(modes)]
        value = (x - 0.25) ** 2 + 0.15 * y + {"A": 0.0, "B": 0.1, "C": 0.3}[mode]
        study.add_trial(
            optuna.trial.create_trial(
                value=value,
                params={"x": x, "y": y, "mode": mode},
                distributions={
                    "x": optuna.distributions.FloatDistribution(-1.0, 1.0),
                    "y": optuna.distributions.FloatDistribution(0.0, 1.0),
                    "mode": optuna.distributions.CategoricalDistribution(modes),
                },
            )
        )
        rows.append({"x": x, "y": y, "mode": mode, "value": value})
    evaluator = optuna.importance.FanovaImportanceEvaluator(
        n_trees=16, max_depth=16, seed=97
    )
    importances = optuna.importance.get_param_importances(
        study, evaluator=evaluator
    )
    mdi = optuna.importance.MeanDecreaseImpurityImportanceEvaluator(
        n_trees=16, max_depth=16, seed=97
    )
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        ped_anova = optuna.importance.PedAnovaImportanceEvaluator(
            target_quantile=0.25, region_quantile=1.0
        )

    def public_result(
        current: optuna.importance.BaseImportanceEvaluator,
        *,
        target: object = None,
    ) -> dict[str, object]:
        values = current.evaluate(study, target=target)
        return {
            "parameter_order": list(values),
            "values": [float(value) for value in values.values()],
        }

    return {
        "trials": rows,
        "evaluator": {"name": "fanova", "n_trees": 16, "max_depth": 16, "seed": 97},
        "parameter_order": list(importances),
        "values": [float(value) for value in importances.values()],
        "public_evaluators": {
            "base_construction_error": "TypeError",
            "fanova": public_result(evaluator),
            "fanova_target": public_result(
                evaluator, target=lambda trial: float(trial.params["y"])
            ),
            "mdi": public_result(mdi),
            "ped_anova": public_result(ped_anova),
        },
    }


def _terminator_contract() -> dict[str, object]:
    minimizing = [5.0, 4.0, 4.5, 4.2]
    maximizing = [1.0, 3.0, 2.5, 2.0]

    def frozen(values: list[float]) -> list[optuna.trial.FrozenTrial]:
        return [optuna.trial.create_trial(value=value) for value in values]

    evaluator = optuna.terminator.BestValueStagnationEvaluator(
        max_stagnation_trials=3
    )
    remaining_minimize = evaluator.evaluate(
        frozen(minimizing), optuna.study.StudyDirection.MINIMIZE
    )
    remaining_maximize = evaluator.evaluate(
        frozen(maximizing), optuna.study.StudyDirection.MAXIMIZE
    )

    callback_study = optuna.create_study()
    callback_study.optimize(
        lambda trial: float(trial.number),
        n_trials=10,
        callbacks=[optuna.study.MaxTrialsCallback(3)],
    )

    stagnation_study = optuna.create_study()
    sequence = iter(minimizing + [9.0, 10.0])
    terminator = optuna.terminator.Terminator(
        improvement_evaluator=optuna.terminator.BestValueStagnationEvaluator(2),
        error_evaluator=optuna.terminator.StaticErrorEvaluator(0.0),
        min_n_trials=1,
    )
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", FutureWarning)
        stagnation_study.optimize(
            lambda _: next(sequence),
            n_trials=6,
            callbacks=[optuna.terminator.TerminatorCallback(terminator)],
        )
    return {
        "remaining_minimize": float(remaining_minimize),
        "remaining_maximize": float(remaining_maximize),
        "max_trials_callback_count": len(callback_study.trials),
        "terminator_trial_count": len(stagnation_study.trials),
        "terminator_values": [float(trial.value) for trial in stagnation_study.trials],
    }


def build_oracle() -> dict[str, object]:
    optuna.logging.set_verbosity(optuna.logging.WARNING)
    if optuna.__version__ != EXPECTED_VERSION:
        raise RuntimeError(
            f"Expected optuna=={EXPECTED_VERSION}, found {optuna.__version__}."
        )
    single = optuna.create_study()
    multi = optuna.create_study(directions=["minimize", "maximize"])
    return {
        "schema": "radia.test.optuna-upstream-oracle.v1",
        "optuna_version": optuna.__version__,
        "numpy_version": np.__version__,
        "scipy_version": scipy.__version__,
        "python_version": sys.version.split()[0],
        "torch_version": importlib.metadata.version("torch"),
        "cmaes_version": importlib.metadata.version("cmaes"),
        "defaults": {
            "anonymous_prefix": "no-name-",
            "single_sampler": type(single.sampler).__name__,
            "multi_sampler": type(multi.sampler).__name__,
            "multi_population_size": multi.sampler._population_size,
            "pruner": type(single.pruner).__name__,
        },
        "exceptions": _exception_contract(),
        "logging": _logging_contract(),
        "sampler_seed_defaults": _sampler_seed_default_contract(),
        "sampler_public_members": _sampler_public_member_contract(),
        "tell": _tell_contract(),
        "trial_pruned_exception": _trial_pruned_exception_contract(),
        "lifecycle_errors": _lifecycle_error_contract(),
        "fixed_trial": _fixed_trial_contract(),
        "frozen_trial": _frozen_trial_contract(),
        "study_management": _study_management_contract(),
        "trials_dataframe": _trials_dataframe_contract(),
        "importance": _importance_contract(),
        "terminator": _terminator_contract(),
        "numpy_random_state_seed_contract": _random_state_contract(),
        "core_api": _core_api_contract(),
        "distributions": _distribution_contract(),
        "distribution_public_members": _distribution_public_member_contract(),
        "search_space": _search_space_contract(),
        "enums": _enum_contract(),
        "unfinished_trials": _unfinished_trial_contract(),
        "base_trial": _base_trial_contract(),
        "base_components": _base_component_contract(),
        "qmc_warnings": _qmc_warning_contract(),
        "distribution_json": _distribution_json_contract(),
        "random_sampler_seed_123": _random_trials(),
        "tpe_sampler_seed_37": _tpe_trials(),
        "custom_tpe_sampler_gamma_weights": _custom_tpe_trials(),
        "tpe_group": _tpe_group_contract(),
        "tpe_categorical_distance": _tpe_categorical_distance_trials(),
        "multiobjective_tpe_sampler_seed_41": _multiobjective_tpe_trials(),
        "mixed_tpe_sampler_seed_43": _mixed_tpe_trials(),
        "grid_sampler_seed_17": _grid_trials(),
        "nsgaii_sampler_seed_19": _nsgaii_trials(),
        "nsgaiii_sampler_seed_23": _nsgaiii_trials(),
        "nsgaii_crossovers_seed_73": _nsgaii_crossover_trials(),
        "brute_force_sampler_seed_29": _brute_force_trials(),
        "conditional_brute_force_sampler_seed_79": _conditional_brute_force_trials(),
        "cmaes_sampler_seed_31": _cmaes_trials(),
        "cmaes_independent_sampler_seed_31": _cmaes_independent_sampler_trials(),
        "scrambled_qmc_sampler_seed_47": _scrambled_qmc_trials(),
        "gp_sampler_seed_53": _gp_trials(),
        "gp_constraints_sampler_seed_89": _gp_constraint_trials(),
        "partial_fixed_sampler_seed_61": _partial_fixed_trials(),
        "multivariate_tpe_sampler_seed_67": _multivariate_tpe_trials(),
        "unscrambled_qmc_sampler_seed_71": _unscrambled_qmc_trials(),
        "pruners": _pruner_contract(),
        "constraints": _constraint_contract(),
    }


def main() -> None:
    destination = Path(__file__).with_name("optuna49_oracle.json")
    destination.write_text(
        json.dumps(build_oracle(), indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(destination)


if __name__ == "__main__":
    main()
