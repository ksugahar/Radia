"""Move a study between MATLAB's table storage and an Optuna storage.

MATLAB keeps its optimization history in tables and a MAT-file, which is not
an ``optuna.storages`` backend, so a ``radia.optuna`` study cannot simply be
opened from Python. That gap is what makes "if MATLAB is missing a feature,
just run Python Optuna" untrue in practice.

This module closes it with an explicit handoff. ``radia.optuna.export_study``
writes a ``radia.optuna.study-export.v1`` document; :func:`into_study` replays
it into any Optuna storage so the Python-only surface (dashboard,
visualization, parameter importance, integration samplers) applies to the same
trials; :func:`from_study` writes the document back for
``radia.optuna.import_study``.

The document carries names as records rather than JSON object keys because a
parameter or attribute name need not be a valid MATLAB field name, and the
handoff must not silently rename anything. Distributions travel as Optuna's
own ``distribution_to_json`` strings, so neither side reinterprets them.
"""

from __future__ import annotations

import argparse
from datetime import datetime
import json
from pathlib import Path
from typing import TYPE_CHECKING, Any
import warnings


if TYPE_CHECKING:  # pragma: no cover - typing only
    import optuna


SCHEMA = "radia.optuna.study-export.v1"
_TIMESTAMP_FORMAT = "%Y-%m-%dT%H:%M:%S.%f"


def _optuna():
    try:
        import optuna  # noqa: PLC0415
    except ModuleNotFoundError as error:  # pragma: no cover - environment
        raise ModuleNotFoundError(
            "The storage bridge needs upstream Optuna. Install it with "
            "'pip install radia-optuna[upstream]'."
        ) from error
    return optuna


def _records(container: dict[str, Any], field: str) -> list[dict[str, Any]]:
    value = container.get(field) or []
    if isinstance(value, dict):  # a single record decoded without its list
        return [value]
    return list(value)


def _parse_timestamp(text: str | None) -> datetime | None:
    if not text:
        return None
    return datetime.strptime(text, _TIMESTAMP_FORMAT)


def _format_timestamp(value: datetime | None) -> str | None:
    if value is None:
        return None
    if value.tzinfo is not None:
        value = value.astimezone().replace(tzinfo=None)
    return value.strftime(_TIMESTAMP_FORMAT)


def load_export(path: str | Path) -> dict[str, Any]:
    """Read and validate an export document."""
    payload = json.loads(Path(path).read_text(encoding="utf-8"))
    if payload.get("schema") != SCHEMA:
        raise ValueError(
            f"Unsupported export schema {payload.get('schema')!r}; expected {SCHEMA!r}."
        )
    return payload


def save_export(payload: dict[str, Any], path: str | Path) -> Path:
    """Write an export document for radia.optuna.import_study."""
    destination = Path(path)
    destination.write_text(
        json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )
    return destination


def frozen_trials(payload: dict[str, Any]) -> list["optuna.trial.FrozenTrial"]:
    """Rebuild every exported trial as an Optuna FrozenTrial."""
    optuna = _optuna()
    trials = []
    for record in sorted(_records(payload, "trials"), key=lambda r: r["number"]):
        params: dict[str, Any] = {}
        distributions: dict[str, Any] = {}
        for item in _records(record, "params"):
            params[item["name"]] = item["value"]
            distributions[item["name"]] = optuna.distributions.json_to_distribution(
                item["distribution"]
            )
        intermediate = {
            int(item["step"]): float(item["value"])
            for item in _records(record, "intermediate_values")
        }
        values = list(record.get("values") or [])
        trial = optuna.trial.create_trial(
            state=optuna.trial.TrialState[record["state"]],
            values=values or None,
            params=params,
            distributions=distributions,
            user_attrs=_attributes(record, "user_attrs"),
            system_attrs=_attributes(record, "system_attrs"),
            intermediate_values=intermediate,
        )
        started = _parse_timestamp(record.get("datetime_start"))
        if started is not None:
            trial.datetime_start = started
        completed = _parse_timestamp(record.get("datetime_complete"))
        if completed is not None:
            trial.datetime_complete = completed
        # Optuna requires a finished trial to carry a completion timestamp,
        # so leave the one create_trial assigned rather than clearing it.
        trials.append(trial)
    return trials


def _attributes(container: dict[str, Any], field: str) -> dict[str, Any]:
    return {
        item["name"]: json.loads(item["value_json"])
        for item in _records(container, field)
    }


def into_study(
    payload: dict[str, Any],
    *,
    storage: Any = None,
    study_name: str | None = None,
    load_if_exists: bool = False,
) -> "optuna.Study":
    """Replay an export into a new Optuna study on the given storage."""
    optuna = _optuna()
    study = optuna.create_study(
        study_name=study_name or payload.get("study_name"),
        storage=storage,
        directions=list(payload["directions"]),
        load_if_exists=load_if_exists,
    )
    for name, value in _attributes(payload, "user_attrs").items():
        study.set_user_attr(name, value)
    metric_names = list(payload.get("metric_names") or [])
    if metric_names:
        # set_metric_names is still flagged experimental upstream; carrying
        # the names across is deliberate, so do not surface that warning to
        # the caller as if the bridge were misusing the API.
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", optuna.exceptions.ExperimentalWarning)
            study.set_metric_names(metric_names)
    study.add_trials(frozen_trials(payload))
    return study


def from_study(study: "optuna.Study") -> dict[str, Any]:
    """Serialize an Optuna study back into an export document."""
    optuna = _optuna()
    records = []
    for trial in study.get_trials(deepcopy=False):
        records.append(
            {
                "number": trial.number,
                "state": trial.state.name,
                "values": list(trial.values) if trial.values is not None else [],
                "params": [
                    {
                        "name": name,
                        "value": value,
                        "distribution": optuna.distributions.distribution_to_json(
                            trial.distributions[name]
                        ),
                    }
                    for name, value in trial.params.items()
                ],
                "user_attrs": _attribute_records(trial.user_attrs),
                "system_attrs": _attribute_records(trial.system_attrs),
                "intermediate_values": [
                    {"step": step, "value": value}
                    for step, value in sorted(trial.intermediate_values.items())
                ],
                "datetime_start": _format_timestamp(trial.datetime_start),
                "datetime_complete": _format_timestamp(trial.datetime_complete),
            }
        )
    return {
        "schema": SCHEMA,
        "study_name": study.study_name,
        "directions": [direction.name.lower() for direction in study.directions],
        "metric_names": list(study.metric_names or []),
        "user_attrs": _attribute_records(study.user_attrs),
        # Study-level system attributes are deprecated upstream since v3.1
        # and slated for removal in v5, so the bridge does not read them.
        "system_attrs": [],
        "trial_count": len(records),
        "trials": records,
    }


def _attribute_records(attributes: dict[str, Any]) -> list[dict[str, Any]]:
    return [
        {"name": name, "value_json": json.dumps(value)}
        for name, value in attributes.items()
    ]


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    subcommands = parser.add_subparsers(dest="command", required=True)

    load = subcommands.add_parser(
        "load", help="replay a MATLAB export into an Optuna storage"
    )
    load.add_argument("export", type=Path)
    load.add_argument("--storage", default=None)
    load.add_argument("--study-name", default=None)
    load.add_argument("--load-if-exists", action="store_true")

    dump = subcommands.add_parser(
        "dump", help="write an Optuna study back out for radia.optuna.import_study"
    )
    dump.add_argument("output", type=Path)
    dump.add_argument("--storage", required=True)
    dump.add_argument("--study-name", required=True)

    args = parser.parse_args()
    optuna = _optuna()
    if args.command == "load":
        study = into_study(
            load_export(args.export),
            storage=args.storage,
            study_name=args.study_name,
            load_if_exists=args.load_if_exists,
        )
        print(
            json.dumps(
                {
                    "schema": "radia-optuna.bridge-load.v1",
                    "study_name": study.study_name,
                    "storage": args.storage,
                    "trials": len(study.get_trials(deepcopy=False)),
                },
                indent=2,
            )
        )
    else:
        study = optuna.load_study(study_name=args.study_name, storage=args.storage)
        destination = save_export(from_study(study), args.output)
        print(
            json.dumps(
                {
                    "schema": "radia-optuna.bridge-dump.v1",
                    "study_name": study.study_name,
                    "output": str(destination),
                    "trials": len(study.get_trials(deepcopy=False)),
                },
                indent=2,
            )
        )


if __name__ == "__main__":  # pragma: no cover
    main()
