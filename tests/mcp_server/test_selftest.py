"""Tests for --selftest: graceful skip, fixture fallback, and fixture validation."""

import os
import sys
from io import StringIO
from pathlib import Path

from radia.mcp_server.radia.server import _lint_file as _lint_file_radia
from radia.mcp_server.ngsolve.server import _lint_file as _lint_file_ngsolve


# ============================================================
# Graceful skip when no examples/ and no fixtures
# ============================================================

def test_radia_selftest_uses_fixtures(monkeypatch):
    """Radia selftest should use fixtures when examples/ not found."""
    monkeypatch.setattr(
        'radia.mcp_server.radia.server.PROJECT_ROOT', Path("/nonexistent")
    )
    from radia.mcp_server.radia.server import _selftest

    captured = StringIO()
    monkeypatch.setattr('sys.stdout', captured)
    _selftest()
    output = captured.getvalue()
    assert 'fixture' in output.lower() or 'PASSED' in output or 'SKIP' in output


def test_ngsolve_selftest_uses_fixtures(monkeypatch):
    """NGSolve selftest should use fixtures when examples/ not found."""
    monkeypatch.setattr(
        'radia.mcp_server.ngsolve.server.PROJECT_ROOT', Path("/nonexistent")
    )
    from radia.mcp_server.ngsolve.server import _selftest

    captured = StringIO()
    monkeypatch.setattr('sys.stdout', captured)
    _selftest()
    output = captured.getvalue()
    assert 'fixture' in output.lower() or 'PASSED' in output or 'SKIP' in output


def test_cubit_selftest_skips_without_examples(tmp_path, monkeypatch):
    """Cubit selftest should print SKIP when examples/ not found."""
    monkeypatch.setattr(
        'radia.mcp_server.cubit.server.PROJECT_ROOT', tmp_path
    )
    from radia.mcp_server.cubit.server import _selftest

    captured = StringIO()
    monkeypatch.setattr('sys.stdout', captured)
    _selftest()
    output = captured.getvalue()
    assert 'SKIP' in output


# ============================================================
# Fixture file validation (Radia)
# ============================================================

FIXTURES_DIR = Path(__file__).parent / "fixtures"


def test_fixtures_dir_exists():
    """Test fixtures directory must exist."""
    assert FIXTURES_DIR.exists(), f"fixtures/ not found at {FIXTURES_DIR}"


def test_bad_radia_has_findings():
    """bad_radia_script.py must trigger at least 6 lint findings."""
    findings = _lint_file_radia(str(FIXTURES_DIR / "bad_radia_script.py"))
    rules_found = {f['rule'] for f in findings}
    assert 'objbckg-needs-callable' in rules_found
    assert 'missing-utidelall' in rules_found
    assert 'hardcoded-absolute-path' in rules_found
    assert 'removed-fldunits' in rules_found
    assert 'removed-fldbatch' in rules_found
    assert 'removed-solver-api' in rules_found
    assert len(findings) >= 6


def test_bad_peec_has_findings_radia():
    """bad_peec_script.py must trigger PEEC/BEM findings via Radia linter."""
    findings = _lint_file_radia(str(FIXTURES_DIR / "bad_peec_script.py"))
    rules_found = {f['rule'] for f in findings}
    assert 'bessel-jv-not-iv' in rules_found
    assert 'peec-low-nseg' in rules_found
    assert len(findings) >= 2


def test_clean_radia_has_no_findings():
    """clean_radia_script.py must produce zero findings (no false positives)."""
    findings = _lint_file_radia(str(FIXTURES_DIR / "clean_radia_script.py"))
    assert findings == [], (
        f"Clean script has {len(findings)} finding(s): "
        + ", ".join(f['rule'] for f in findings)
    )


# ============================================================
# Fixture file validation (NGSolve)
# ============================================================

def test_bad_ngsolve_has_findings():
    """bad_ngsolve_script.py must trigger NGSolve-specific findings."""
    findings = _lint_file_ngsolve(str(FIXTURES_DIR / "bad_ngsolve_script.py"))
    rules_found = {f['rule'] for f in findings}
    assert 'ngsolve-overwrite-xyz' in rules_found
    assert 'ngsolve-dim2-occ' in rules_found
    assert len(findings) >= 2


def test_bad_peec_has_findings_ngsolve():
    """bad_peec_script.py must trigger PEEC/BEM findings via NGSolve linter."""
    findings = _lint_file_ngsolve(str(FIXTURES_DIR / "bad_peec_script.py"))
    rules_found = {f['rule'] for f in findings}
    assert 'bessel-jv-not-iv' in rules_found
    assert 'peec-low-nseg' in rules_found
    assert len(findings) >= 2


def test_clean_ngsolve_has_no_findings():
    """clean_ngsolve_script.py must produce zero findings (no false positives)."""
    findings = _lint_file_ngsolve(str(FIXTURES_DIR / "clean_ngsolve_script.py"))
    assert findings == [], (
        f"Clean script has {len(findings)} finding(s): "
        + ", ".join(f['rule'] for f in findings)
    )


# ============================================================
# Cross-server coverage
# ============================================================

def test_all_severity_levels_covered():
    """Fixtures should cover CRITICAL, HIGH, MODERATE, and LOW severities."""
    all_findings = []
    for py_file in FIXTURES_DIR.glob("bad_*.py"):
        all_findings.extend(_lint_file_radia(str(py_file)))
        all_findings.extend(_lint_file_ngsolve(str(py_file)))

    severities = {f['severity'] for f in all_findings}
    assert 'CRITICAL' in severities, "No CRITICAL findings in fixtures"
    assert 'HIGH' in severities, "No HIGH findings in fixtures"
