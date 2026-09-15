"""Cubit selftest and repository audit contracts."""
import contextlib
import errno
import io


def test_cubit_selftest_skips_repo_audit_by_default(monkeypatch, tmp_path):
    from cubit_mesh_export.mcp import server

    (tmp_path / "docs").mkdir()
    monkeypatch.setattr(server, "PROJECT_ROOT", tmp_path)

    out = io.StringIO()
    with contextlib.redirect_stdout(out):
        server._selftest()

    text = out.getvalue()
    assert "repo audit: SKIPPED" in text
    assert "PASSED" in text
    assert "Cubit Export Lint Report" not in text


def test_mesh_selftest_cli_tolerates_closed_stdout(monkeypatch):
    from cubit_mesh_export.mcp import server as cubit_server

    for module in (cubit_server,):
        def raise_closed_pipe(*, audit_repo=False):
            raise OSError(errno.EINVAL, "Invalid argument")

        monkeypatch.setattr(module, "_selftest", raise_closed_pipe)
        monkeypatch.setattr(module.sys, "argv", ["cmd", "--selftest", "--audit-repo"])
        module.main()


def test_mesh_status_tools_expose_selftest_and_audit_commands():
    from cubit_mesh_export.mcp import server as cubit_server

    cubit = cubit_server.mcp._tool_manager._tools["cubit_status"].fn()

    assert cubit["selftest_command"] == "mcp-server-cubit --selftest"
    assert cubit["audit_command"] == "mcp-server-cubit --selftest --audit-repo"
    assert "cubit_status" in cubit["tools"]


def test_mesh_audit_summary_tools_are_machine_readable(monkeypatch, tmp_path):
    from cubit_mesh_export.mcp import server as cubit_server

    examples = tmp_path / "examples"
    examples.mkdir()
    (examples / "a.py").write_text("print('a')\n", encoding="utf-8")
    (examples / "b.py").write_text("print('b')\n", encoding="utf-8")

    def fake_lint(filepath: str):
        if filepath.endswith("a.py"):
            return [
                {"line": 1, "severity": "HIGH", "rule": "alpha", "message": "x"},
                {"line": 2, "severity": "LOW", "rule": "alpha", "message": "y"},
                {"line": 3, "severity": "CRITICAL", "rule": "beta", "message": "z"},
            ]
        return []

    for module, tool_name in (
        (cubit_server, "cubit_audit_summary"),
    ):
        monkeypatch.setattr(module, "PROJECT_ROOT", tmp_path)
        monkeypatch.setattr(module, "_lint_file", fake_lint)
        summary = getattr(module, tool_name)("examples", top_n=1)
        assert summary["ok"] is True
        assert summary["files_scanned"] == 2
        assert summary["files_with_findings"] == 1
        assert summary["total_findings"] == 3
        assert summary["clean"] is False
        assert summary["by_severity"]["HIGH"] == 1
        assert summary["by_severity"]["LOW"] == 1
        assert summary["by_severity"]["CRITICAL"] == 1
        assert summary["top_rules"] == [{
            "rule": "alpha",
            "count": 2,
            "action": "Inspect representative findings and add a specific remediation note.",
        }]
        assert summary["dominant_rule"] == summary["top_rules"][0]
        assert summary["top_files"][0]["path"] == "examples/a.py"
