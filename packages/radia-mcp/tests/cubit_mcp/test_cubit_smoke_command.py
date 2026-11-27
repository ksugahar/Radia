"""Guards: the export smoke test drives the isolated export worker on the bundled
sample and fails loudly on a runtime or export failure before any vol_check.

Assumes: no Cubit, Netgen or worker run; ``export_runner.run_export`` is
replaced by a recorder.  The solver-ready gate itself is radia.vol_check's.
"""
from __future__ import annotations

from pathlib import Path

import pytest

from radia_mcp.cubit import export_runner, smoke_test


def test_bundled_sample_is_the_default_journal():
    sample = smoke_test._find_sample_jou()
    assert sample.name == "solver_ready_sample.jou" and sample.parent.name == "gui", sample
    assert sample.is_file()
    with pytest.raises(RuntimeError, match="does not exist"):
        smoke_test._find_sample_jou(str(Path(sample).with_name("missing.jou")))


def test_smoke_exports_through_the_worker_with_the_given_runtime(monkeypatch, tmp_path, capsys):
    calls = []

    def run_export(journal, out_path, order, **kwargs):
        calls.append((Path(journal).name, Path(out_path).name, order, kwargs))
        return {"status": "error", "error_type": "RuntimeError", "error": "boom",
                "stderr_tail": "worker said no"}

    monkeypatch.setattr(export_runner, "run_export", run_export)
    monkeypatch.setattr(smoke_test.tempfile, "mkdtemp", lambda **k: str(tmp_path / "work"))
    (tmp_path / "work").mkdir()
    code = smoke_test.run_smoke_test(order=3, plugin_dir="P", cubit_bin="B", python="Y", timeout=7)
    assert code == 1
    assert calls == [("solver_ready_sample.jou", "smoke.vol", 3,
                      {"cubit_bin": "B", "python": "Y", "plugin_dir": "P", "timeout_s": 7})], calls
    out = capsys.readouterr().out
    assert "[FAIL] export failed: RuntimeError: boom" in out and "worker said no" in out, out


def test_missing_runtime_is_a_loud_failure(monkeypatch, tmp_path, capsys):
    def run_export(*args, **kwargs):
        raise export_runner.WorkerRuntimeError("Native plugin directory not set")

    monkeypatch.setattr(export_runner, "run_export", run_export)
    monkeypatch.setattr(smoke_test.tempfile, "mkdtemp", lambda **k: str(tmp_path / "work"))
    (tmp_path / "work").mkdir()
    assert smoke_test.run_smoke_test() == 1
    assert "[FAIL] export runtime: WorkerRuntimeError: Native plugin directory not set" in (
        capsys.readouterr().out)


def test_solver_ready_gate_is_radia_vol_check(monkeypatch, tmp_path):
    import json
    import sys
    import types

    seen = {}
    fake = types.ModuleType("radia.vol_check")
    fake.REPORT_SCHEMA = "cubit-mesh-export.vol-check.v1"

    def check_consistency(path, **kwargs):
        seen.update(path=Path(path), **kwargs)
        return {"schema": fake.REPORT_SCHEMA, "passed": False, "warnings": ["fake verdict"]}

    fake.check_consistency = check_consistency
    monkeypatch.setitem(sys.modules, "radia.vol_check", fake)
    vol = tmp_path / "smoke.vol"
    vol.write_text('materials\n1\n1 "air"\n\nbcnames\n1\n1 "outer"\n', encoding="utf-8")
    Path(str(vol) + ".json").write_text(json.dumps(
        {"materials": {}, "boundaries": {}, "edges": {}, "n_elements": 1, "n_points": 4,
         "order": 2, "export_time_s": 0.1}), encoding="utf-8")
    result = smoke_test._validate_exported_vol(vol, order=2, expect=["outer"],
                                               expect_materials=["air"], threshold=1.0)
    assert seen["path"] == vol and seen["strict_labels"] is True, seen
    assert seen["required_boundaries"] == ("outer",) and seen["required_materials"] == ("air",), seen
    assert seen["min_curve_order"] == 2 and seen["require_tetrahedra"] is True, seen
    assert result["passed"] is False and "fake verdict" in result["issues"], result
    assert json.loads((tmp_path / "vol-check.json").read_text(encoding="utf-8"))["passed"] is False
