"""Keep lightweight MCP selftests separate from repo-wide audits."""

from __future__ import annotations

import contextlib
import errno
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


def test_mesh_selftest_cli_tolerates_closed_stdout(monkeypatch):
    from radia_mcp.cubit import server as cubit_server
    from radia_mcp.gmsh import server as gmsh_server

    for module in (cubit_server, gmsh_server):
        def raise_closed_pipe(*, audit_examples=False):
            raise OSError(errno.EINVAL, "Invalid argument")

        monkeypatch.setattr(module, "_selftest", raise_closed_pipe)
        monkeypatch.setattr(module.sys, "argv", ["cmd", "--selftest", "--audit-examples"])
        module.main()
