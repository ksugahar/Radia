from __future__ import annotations

import json
from pathlib import Path
import sys

import pytest

PACKAGE_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PACKAGE_ROOT / "src"))

from radia_optuna import bridge  # noqa: E402


def _payload() -> dict:
    distribution = json.dumps(
        {
            "name": "FloatDistribution",
            "attributes": {"step": None, "low": -1.0, "high": 1.0, "log": False},
        }
    )
    return {
        "schema": bridge.SCHEMA,
        "study_name": "bridge-demo",
        "directions": ["minimize"],
        "metric_names": ["loss"],
        "user_attrs": [{"name": "owner-name", "value_json": '"radia"'}],
        "system_attrs": [{"name": "source-id", "value_json": "7"}],
        "trial_count": 2,
        "trials": [
            {
                "number": 0,
                "state": "COMPLETE",
                "values": [0.25],
                "params": [
                    {"name": "x-1", "value": 0.5, "distribution": distribution}
                ],
                "user_attrs": [{"name": "trial-tag", "value_json": '"t1"'}],
                "system_attrs": [{"name": "worker-id", "value_json": "3"}],
                "intermediate_values": [{"step": 0, "value": 1.5}],
                "constraint_present": True,
                "constraints": [-0.5, 0.0],
                "datetime_start": "2026-08-28T10:00:00.000000",
                "datetime_complete": "2026-08-28T10:00:01.500000",
            },
            {
                "number": 1,
                "state": "FAIL",
                "values": [],
                "params": [],
                "user_attrs": [],
                "system_attrs": [],
                "intermediate_values": [],
                "constraint_present": False,
                "constraints": [],
                "datetime_start": "2026-08-28T10:00:02.000000",
                "datetime_complete": "2026-08-28T10:00:02.100000",
            },
        ],
    }


def test_schema_validation(tmp_path):
    path = tmp_path / "wrong.json"
    path.write_text(json.dumps({"schema": "wrong"}), encoding="utf-8")
    with pytest.raises(ValueError, match="Unsupported export schema"):
        bridge.load_export(path)


def test_roundtrip_through_real_sqlite_storage(tmp_path):
    optuna = pytest.importorskip("optuna")
    optuna.logging.set_verbosity(optuna.logging.WARNING)
    source = _payload()
    storage = f"sqlite:///{(tmp_path / 'bridge.db').as_posix()}"
    study = bridge.into_study(source, storage=storage)
    reopened = optuna.load_study(study_name=study.study_name, storage=storage)
    echoed = bridge.from_study(reopened)

    assert echoed["directions"] == source["directions"]
    assert echoed["metric_names"] == ["loss"]
    assert echoed["user_attrs"] == source["user_attrs"]
    assert echoed["system_attrs"] == source["system_attrs"]
    first = echoed["trials"][0]
    assert first["params"][0]["name"] == "x-1"
    assert first["user_attrs"] == source["trials"][0]["user_attrs"]
    assert first["system_attrs"] == source["trials"][0]["system_attrs"]
    assert first["constraint_present"] is True
    assert first["constraints"] == [-0.5, 0.0]
    assert first["intermediate_values"] == [{"step": 0, "value": 1.5}]
    assert first["datetime_start"] == "2026-08-28T10:00:00.000000"
    assert first["datetime_complete"] == "2026-08-28T10:00:01.500000"


def test_bridge_refuses_an_unpinned_optuna(monkeypatch):
    import optuna

    monkeypatch.setattr(optuna, "__version__", "4.8.0")
    with pytest.raises(RuntimeError, match="requires optuna==4.9.0"):
        bridge._optuna()
