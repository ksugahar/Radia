"""The MATLAB table storage is not an optuna.storages backend.

Without a handoff, "if MATLAB is missing a feature, run Python Optuna" is not
actually true: the study cannot be opened from Python at all. These tests lock
the handoff document instead, so an export produced by
``radia.optuna.export_study`` replays into a real Optuna storage and comes back
unchanged.
"""

from __future__ import annotations

import json
from pathlib import Path
import sys

import pytest


PACKAGE_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PACKAGE_ROOT / "src"))

from radia_optuna import bridge  # noqa: E402


def _export() -> dict:
    """An export shaped exactly like radia.optuna.export_study writes one."""
    float_json = json.dumps(
        {
            "name": "FloatDistribution",
            "attributes": {"step": None, "low": -1.0, "high": 1.0, "log": False},
        }
    )
    categorical_json = json.dumps(
        {"name": "CategoricalDistribution", "attributes": {"choices": ["a", "b"]}}
    )
    return {
        "schema": bridge.SCHEMA,
        "study_name": "bridge-demo",
        "directions": ["minimize", "minimize"],
        "metric_names": ["loss", "cost"],
        "user_attrs": [{"name": "owner", "value_json": '"radia"'}],
        "system_attrs": [],
        "trial_count": 3,
        "trials": [
            {
                "number": 0,
                "state": "COMPLETE",
                "values": [0.25, 1.0],
                # "x-1" is not a valid MATLAB field name; the document must
                # carry the real name, not an escaped one.
                "params": [
                    {"name": "x-1", "value": 0.5, "distribution": float_json},
                    {"name": "mode", "value": "a", "distribution": categorical_json},
                ],
                "user_attrs": [{"name": "tag", "value_json": '"t1"'}],
                "system_attrs": [],
                "intermediate_values": [
                    {"step": 0, "value": 1.5},
                    {"step": 1, "value": 0.8},
                ],
                "datetime_start": "2026-08-28T10:00:00.000000",
                "datetime_complete": "2026-08-28T10:00:01.500000",
            },
            {
                "number": 1,
                "state": "PRUNED",
                "values": [],
                "params": [
                    {"name": "x-1", "value": -0.25, "distribution": float_json},
                    {"name": "mode", "value": "b", "distribution": categorical_json},
                ],
                "user_attrs": [],
                "system_attrs": [],
                "intermediate_values": [{"step": 0, "value": 9.0}],
                "datetime_start": "2026-08-28T10:00:02.000000",
                "datetime_complete": "2026-08-28T10:00:02.250000",
            },
            {
                "number": 2,
                "state": "FAIL",
                "values": [],
                "params": [
                    {"name": "x-1", "value": 0.75, "distribution": float_json}
                ],
                "user_attrs": [],
                "system_attrs": [],
                "intermediate_values": [],
                "datetime_start": "2026-08-28T10:00:03.000000",
                "datetime_complete": "2026-08-28T10:00:03.100000",
            },
        ],
    }


def test_export_schema_is_rejected_when_it_does_not_match(tmp_path):
    path = tmp_path / "wrong.json"
    path.write_text(json.dumps({"schema": "something.else"}), encoding="utf-8")
    with pytest.raises(ValueError, match="Unsupported export schema"):
        bridge.load_export(path)


def test_export_replays_into_an_optuna_storage_and_comes_back_unchanged(tmp_path):
    optuna = pytest.importorskip("optuna")
    optuna.logging.set_verbosity(optuna.logging.WARNING)
    source = _export()

    storage = f"sqlite:///{(tmp_path / 'bridge.db').as_posix()}"
    study = bridge.into_study(source, storage=storage)

    trials = study.get_trials(deepcopy=False)
    assert [trial.state.name for trial in trials] == ["COMPLETE", "PRUNED", "FAIL"]
    assert trials[0].params == {"x-1": 0.5, "mode": "a"}
    assert trials[0].values == [0.25, 1.0]
    assert trials[0].intermediate_values == {0: 1.5, 1: 0.8}
    assert trials[0].user_attrs == {"tag": "t1"}
    assert study.user_attrs == {"owner": "radia"}
    assert study.metric_names == ["loss", "cost"]
    assert isinstance(
        trials[0].distributions["mode"], optuna.distributions.CategoricalDistribution
    )

    # The whole point: Python-only surfaces now apply to these trials.
    assert len(study.best_trials) >= 1

    reopened = optuna.load_study(study_name=study.study_name, storage=storage)
    returned = bridge.from_study(reopened)
    assert returned["schema"] == bridge.SCHEMA
    assert returned["directions"] == source["directions"]
    assert returned["trial_count"] == source["trial_count"]

    for original, echoed in zip(source["trials"], returned["trials"]):
        assert echoed["state"] == original["state"]
        assert echoed["values"] == original["values"]
        assert {p["name"]: p["value"] for p in echoed["params"]} == {
            p["name"]: p["value"] for p in original["params"]
        }
        assert echoed["intermediate_values"] == original["intermediate_values"]
        assert echoed["datetime_start"] == original["datetime_start"]
        assert echoed["datetime_complete"] == original["datetime_complete"]
        assert {a["name"] for a in echoed["user_attrs"]} == {
            a["name"] for a in original["user_attrs"]
        }


def test_distributions_travel_as_upstream_json(tmp_path):
    optuna = pytest.importorskip("optuna")
    optuna.logging.set_verbosity(optuna.logging.WARNING)
    study = bridge.into_study(_export())
    trial = study.get_trials(deepcopy=False)[0]
    # Neither side reinterprets the distribution: what MATLAB wrote is exactly
    # what optuna.distributions.json_to_distribution consumed.
    echoed = bridge.from_study(study)["trials"][0]["params"]
    by_name = {item["name"]: item["distribution"] for item in echoed}
    assert json.loads(by_name["x-1"]) == {
        "name": "FloatDistribution",
        "attributes": {"step": None, "low": -1.0, "high": 1.0, "log": False},
    }
    assert (
        optuna.distributions.json_to_distribution(by_name["mode"])
        == trial.distributions["mode"]
    )


def test_a_finished_trial_without_a_completion_time_is_still_accepted():
    """Optuna requires datetime_complete on a finished trial.

    A MATLAB study always supplies one (addTrial falls back to the start
    time), but a hand-written or truncated document may not. The bridge must
    keep the timestamp create_trial assigned instead of clearing it and
    tripping Optuna's own validation.
    """
    optuna = pytest.importorskip("optuna")
    optuna.logging.set_verbosity(optuna.logging.WARNING)
    payload = _export()
    payload["trials"] = [payload["trials"][2]]
    payload["trial_count"] = 1
    payload["trials"][0]["datetime_complete"] = None

    study = bridge.into_study(payload)
    trial = study.get_trials(deepcopy=False)[0]
    assert trial.state is optuna.trial.TrialState.FAIL
    assert trial.datetime_complete is not None
