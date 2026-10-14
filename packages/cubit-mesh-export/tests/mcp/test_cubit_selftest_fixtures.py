"""Exporter-owned selftest and lint fixture regression contracts."""

from io import StringIO
from pathlib import Path

from cubit_mesh_export.mcp.server import _lint_file, _selftest

FIXTURES_DIR = Path(__file__).parent / "fixtures"


def test_cubit_selftest_runs_without_examples(tmp_path, monkeypatch):
    monkeypatch.setattr("cubit_mesh_export.mcp.server.PROJECT_ROOT", tmp_path)
    captured = StringIO()
    monkeypatch.setattr("sys.stdout", captured)
    _selftest()
    output = captured.getvalue()
    assert "PASSED" in output or "SKIP" in output


def test_bad_cubit_has_findings():
    findings = _lint_file(str(FIXTURES_DIR / "bad_cubit_script.py"))
    rules_found = {f["rule"] for f in findings}
    assert "deleted-api-usage" in rules_found
    assert "hardcoded-absolute-path" in rules_found
    assert len(findings) >= 4


def test_clean_cubit_has_no_findings():
    findings = _lint_file(str(FIXTURES_DIR / "clean_cubit_script.py"))
    assert findings == [], (
        f"Clean script has {len(findings)} finding(s): "
        + ", ".join(f["rule"] for f in findings)
    )
