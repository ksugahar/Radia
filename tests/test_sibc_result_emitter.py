"""Artifact regressions without running numerical validation sweeps."""

import importlib.util
import ast
import json
import math
from pathlib import Path
from types import SimpleNamespace

import numpy as np
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


def test_unvalidated_cylinder_candidates_are_pending_not_executed(monkeypatch):
    observed = []
    def load(relative):
        observed.append(relative)
        label = dict(emitter.CASES)[relative]
        return SimpleNamespace(summary=lambda: {"case": label})
    monkeypatch.setattr(emitter, "load", load)
    result = emitter.build()
    for label in ("cylinder_bulk_tower", "cylinder_two_point_ladder"):
        pending = result["pending_cases"][label]
        assert pending["status"] == "HOLD"
        assert pending["script"] not in observed
        assert label not in result["cases"]
        assert "known_result" not in pending


def _measurement_only():
    # Extract only the measurement function, never the script's symbolic or
    # native solver initialization. Numerical values below are test doubles.
    path = PATH.parent / "cylinder/03_rank_N_bulk_sweep.py"
    tree = ast.parse(path.read_text(encoding="utf-8"))
    function = next(node for node in tree.body
                    if isinstance(node, ast.FunctionDef) and node.name == "summary")
    namespace = {"np": np, "math": math, "SWEEP": np.array([1e4, 1e5, 1e6]),
                 "WALL_BAND_HZ": (1e4, 1e6), "A_NUM": 0.005,
                 "SIGMA": 1, "MU": 1, "METRIC": "complex-relative-error"}
    exec(compile(ast.Module(body=[function], type_ignores=[]), str(path), "exec"), namespace)
    return namespace


def test_bulk_summary_preserves_phase_error_rank_and_open_wall_band():
    namespace = _measurement_only()
    namespace["Y_exact_cylinder"] = lambda *args: 1 + 0j
    namespace["Y_mixed_galerkin"] = lambda s, rank: (
        1j if 1e4 < s.imag / (2 * math.pi) < 1e6 else 10 + 0j)
    result = namespace["summary"]()
    assert result["validation_status"] == "HOLD"
    assert result["sweep"]["wall_band_endpoints"] == "excluded"
    for rank, row in result["by_n_bulk"].items():
        assert row["n_unknowns"] == int(rank) + 1
        assert row["wall_band_max_error_pct"] == pytest.approx(100 * math.sqrt(2))
        assert row["max_error_pct"] == pytest.approx(900)


@pytest.mark.parametrize("reference,computed", [(0, 1), (float("nan"), 1), (1, float("inf"))])
def test_bulk_summary_rejects_invalid_admittance(reference, computed):
    namespace = _measurement_only()
    namespace["Y_exact_cylinder"] = lambda *args: reference
    namespace["Y_mixed_galerkin"] = lambda *args: computed
    with pytest.raises(ValueError, match="admittance|reference"):
        namespace["summary"]()


def test_bulk_summary_rejects_empty_wall_band():
    namespace = _measurement_only()
    namespace["SWEEP"] = np.array([1, 10])
    with pytest.raises(ValueError, match="wall band"):
        namespace["summary"]()
