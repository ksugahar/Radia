"""Artifact regressions without running numerical validation sweeps."""

import importlib.util
import json
from pathlib import Path

import pytest


PATH = Path(__file__).resolve().parents[1] / "validation_test/mixed_galerkin/emit_results.py"
SPEC = importlib.util.spec_from_file_location("sibc_result_emitter", PATH)
emitter = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(emitter)


@pytest.mark.parametrize("key", ["max_error_pct", "gamma1_max_error_pct", "planar_max_error_pct"])
def test_zero_error_is_not_replaced(key):
    assert emitter.headline_error({key: 0.0, "by_n_dof": {"4": {"max_error_pct": 5}}}) == 0.0


def test_missing_metric_falls_back():
    assert emitter.headline_error({"max_error_pct": None, "planar_max_error_pct": 2}) == 2
    assert emitter.headline_error({}) is None


@pytest.mark.parametrize("value", [float("nan"), float("inf"), -float("inf")])
def test_nonfinite_result_preserves_existing_artifact(monkeypatch, tmp_path, value):
    output = tmp_path / "result.json"
    output.write_text('{"previous": true}\n', encoding="utf-8")
    monkeypatch.setattr(emitter, "build", lambda: {"cases": {"case": {"max_error_pct": value}}})
    monkeypatch.setattr("sys.argv", ["emit_results.py", "--out", str(output)])
    with pytest.raises(ValueError, match="Out of range float"):
        emitter.main()
    assert json.loads(output.read_text(encoding="utf-8")) == {"previous": True}


def test_valid_result_is_saved_and_zero_is_printed(monkeypatch, tmp_path, capsys):
    output = tmp_path / "nested" / "result.json"
    data = {"cases": {"case": {"max_error_pct": 0.0}}}
    monkeypatch.setattr(emitter, "build", lambda: data)
    monkeypatch.setattr("sys.argv", ["emit_results.py", "--out", str(output)])
    assert emitter.main() == 0
    assert json.loads(output.read_text(encoding="utf-8")) == data
    assert "0.00000 %" in capsys.readouterr().out
