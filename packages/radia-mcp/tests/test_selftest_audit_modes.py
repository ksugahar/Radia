"""Keep lightweight MCP selftests separate from repo-wide audits."""

from __future__ import annotations

import contextlib
import io


def test_cubit_selftest_skips_examples_audit_by_default(monkeypatch, tmp_path):
    from radia_mcp.cubit import server

    (tmp_path / "examples").mkdir()
    monkeypatch.setattr(server, "PROJECT_ROOT", tmp_path)

    out = io.StringIO()
    with contextlib.redirect_stdout(out):
        server._selftest()

    text = out.getvalue()
    assert "examples audit: SKIPPED" in text
    assert "PASSED" in text
    assert "Cubit Export Lint Report" not in text


def test_gmsh_selftest_skips_examples_audit_by_default(monkeypatch, tmp_path):
    from radia_mcp.gmsh import server

    (tmp_path / "examples").mkdir()
    monkeypatch.setattr(server, "PROJECT_ROOT", tmp_path)

    out = io.StringIO()
    with contextlib.redirect_stdout(out):
        server._selftest()

    text = out.getvalue()
    assert "examples audit: SKIPPED" in text
    assert "PASSED" in text
    assert "GMSH issues:" not in text
