"""Explicit handoff between MATLAB study tables and Optuna storages."""

from __future__ import annotations

import argparse
from datetime import datetime
import json
from pathlib import Path
from typing import TYPE_CHECKING, Any
import warnings

if TYPE_CHECKING:  # pragma: no cover
    import optuna

SCHEMA = "radia.optuna.study-export.v1"
_TIMESTAMP_FORMAT = "%Y-%m-%dT%H:%M:%S.%f"
_CONSTRAINTS_KEY = "constraints"


def _optuna():
    try:
        import optuna  # noqa: PLC0415
    except ModuleNotFoundError as error:  # pragma: no cover
        raise ModuleNotFoundError(
            "The storage bridge needs upstream Optuna. Install "
            "'radia-optuna[upstream]'."
        ) from error
    if optuna.__version__ != "4.9.0":
        raise RuntimeError(
            f"The bridge requires optuna==4.9.0, found {optuna.__version__}."
        )
    return optuna


def _records(container: dict[str, Any], field: str) -> list[dict[str, Any]]:
    value = container.get(field) or []
    return [value] if isinstance(value, dict) else list(value)


def _attributes(container: dict[str, Any], field: str) -> dict[str, Any]:
    return {
        item["name"]: json.loads(item["value_json"])
        for item in _records(container, field)
    }


def _attribute_records(attributes: dict[str, Any]) -> list[dict[str, Any]]:
    return [
        {"name": name, "value_json": json.dumps(value, ensure_ascii=False)}
        for name, value in attributes.items()
    ]


def _parse_timestamp(text: str | None) -> datetime | None:
    return datetime.strptime(text, _TIMESTAMP_FORMAT) if text else None


def _format_timestamp(value: datetime | None) -> str | None:
    if value is None:
        return None
    if value.tzinfo is not None:
        value = value.astimezone().replace(tzinfo=None)
    return value.strftime(_TIMESTAMP_FORMAT)


def load_export(path: str | Path) -> dict[str, Any]:
    """Read and validate one handoff document."""
    payload = json.loads(Path(path).read_text(encoding="utf-8"))
    if payload.get("schema") != SCHEMA:
        raise ValueError(
            f"Unsupported export schema {payload.get('schema')!r}; "
            f"expected {SCHEMA!r}."
        )
    return payload


def save_export(payload: dict[str, Any], path: str | Path) -> Path:
    """Write deterministic UTF-8 JSON for radia.optuna.import_study."""
    destination = Path(path)
    destination.write_text(
        json.dumps(payload, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
    return destination


def frozen_trials(payload: dict[str, Any]) -> list["optuna.trial.FrozenTrial"]:
    """Rebuild exported records as upstream FrozenTrial objects."""
    optuna = _optuna()
    trials = []
    for record in sorted(_records(payload, "trials"), key=lambda item: item["number"]):
        params: dict[str, Any] = {}
        distributions: dict[str, Any] = {}
        for item in _records(record, "params"):
            params[item["name"]] = item["value"]
            distributions[item["name"]] = (
                optuna.distributions.json_to_distribution(item["distribution"])
            )
        system_attrs = _attributes(record, "system_attrs")
        if record.get("constraint_present", False):
            system_attrs[_CONSTRAINTS_KEY] = list(record.get("constraints") or [])
        trial = optuna.trial.create_trial(
            state=optuna.trial.TrialState[record["state"]],
            values=list(record.get("values") or []) or None,
            params=params,
            distributions=distributions,
            user_attrs=_attributes(record, "user_attrs"),
            system_attrs=system_attrs,
            intermediate_values={
                int(item["step"]): float(item["value"])
                for item in _records(record, "intermediate_values")
            },
        )
        started = _parse_timestamp(record.get("datetime_start"))
        if started is not None:
            trial.datetime_start = started
        completed = _parse_timestamp(record.get("datetime_complete"))
        if completed is not None:
            trial.datetime_complete = completed
        trials.append(trial)
    return trials


def into_study(
    payload: dict[str, Any],
    *,
    storage: Any = None,
    study_name: str | None = None,
    load_if_exists: bool = False,
) -> "optuna.Study":
    """Replay an export into an upstream Optuna storage."""
    optuna = _optuna()
    study = optuna.create_study(
        study_name=study_name or payload.get("study_name"),
        storage=storage,
        directions=list(payload["directions"]),
        load_if_exists=load_if_exists,
    )
    for name, value in _attributes(payload, "user_attrs").items():
        study.set_user_attr(name, value)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", FutureWarning)
        for name, value in _attributes(payload, "system_attrs").items():
            study.set_system_attr(name, value)
    metric_names = list(payload.get("metric_names") or [])
    if metric_names:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", optuna.exceptions.ExperimentalWarning)
            study.set_metric_names(metric_names)
    study.add_trials(frozen_trials(payload))
    return study


def from_study(study: "optuna.Study") -> dict[str, Any]:
    """Serialize an upstream study back into a MATLAB handoff document."""
    optuna = _optuna()
    records = []
    for trial in study.get_trials(deepcopy=False):
        system_attrs = dict(trial.system_attrs)
        constraint_present = _CONSTRAINTS_KEY in system_attrs
        constraints = list(system_attrs.pop(_CONSTRAINTS_KEY, []))
        records.append(
            {
                "number": trial.number,
                "state": trial.state.name,
                "values": list(trial.values) if trial.values is not None else [],
                "params": [
                    {
                        "name": name,
                        "value": value,
                        "distribution": (
                            optuna.distributions.distribution_to_json(
                                trial.distributions[name]
                            )
                        ),
                    }
                    for name, value in trial.params.items()
                ],
                "user_attrs": _attribute_records(trial.user_attrs),
                "system_attrs": _attribute_records(system_attrs),
                "intermediate_values": [
                    {"step": step, "value": value}
                    for step, value in sorted(trial.intermediate_values.items())
                ],
                "constraint_present": constraint_present,
                "constraints": constraints,
                "datetime_start": _format_timestamp(trial.datetime_start),
                "datetime_complete": _format_timestamp(trial.datetime_complete),
            }
        )
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", FutureWarning)
        system_attrs = dict(study.system_attrs)
    # Optuna persists metric names in this private system-attribute key;
    # metric_names already carries the public value in the handoff schema.
    system_attrs.pop("study:metric_names", None)
    return {
        "schema": SCHEMA,
        "study_name": study.study_name,
        "directions": [direction.name.lower() for direction in study.directions],
        "metric_names": list(study.metric_names or []),
        "user_attrs": _attribute_records(study.user_attrs),
        "system_attrs": _attribute_records(system_attrs),
        "trial_count": len(records),
        "trials": records,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    subcommands = parser.add_subparsers(dest="command", required=True)
    load = subcommands.add_parser("load", help="load a MATLAB export")
    load.add_argument("export", type=Path)
    load.add_argument("--storage", default=None)
    load.add_argument("--study-name", default=None)
    load.add_argument("--load-if-exists", action="store_true")
    dump = subcommands.add_parser("dump", help="dump an Optuna study")
    dump.add_argument("output", type=Path)
    dump.add_argument("--storage", required=True)
    dump.add_argument("--study-name", required=True)
    args = parser.parse_args()
    if args.command == "load":
        study = into_study(
            load_export(args.export),
            storage=args.storage,
            study_name=args.study_name,
            load_if_exists=args.load_if_exists,
        )
        result = {
            "schema": "radia-optuna.bridge-load.v1",
            "study_name": study.study_name,
            "trials": len(study.get_trials(deepcopy=False)),
        }
    else:
        optuna = _optuna()
        study = optuna.load_study(study_name=args.study_name, storage=args.storage)
        destination = save_export(from_study(study), args.output)
        result = {
            "schema": "radia-optuna.bridge-dump.v1",
            "study_name": study.study_name,
            "output": str(destination),
            "trials": len(study.get_trials(deepcopy=False)),
        }
    print(json.dumps(result, indent=2))


if __name__ == "__main__":  # pragma: no cover
    main()
