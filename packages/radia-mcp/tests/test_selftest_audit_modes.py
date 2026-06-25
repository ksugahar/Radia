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


def test_mesh_status_tools_expose_selftest_and_audit_commands():
    from radia_mcp.cubit import server as cubit_server
    from radia_mcp.gmsh import server as gmsh_server

    cubit = cubit_server.mcp._tool_manager._tools["cubit_status"].fn()
    gmsh = gmsh_server.mcp._tool_manager._tools["gmsh_status"].fn()

    assert cubit["selftest_command"] == "mcp-server-cubit --selftest"
    assert cubit["audit_command"] == "mcp-server-cubit --selftest --audit-examples"
    assert "cubit_status" in cubit["tools"]
    assert gmsh["selftest_command"] == "mcp-server-gmsh --selftest"
    assert gmsh["audit_command"] == "mcp-server-gmsh --selftest --audit-examples"
    assert "gmsh_status" in gmsh["tools"]


def test_mesh_audit_summary_tools_are_machine_readable(monkeypatch, tmp_path):
    from radia_mcp.cubit import server as cubit_server
    from radia_mcp.gmsh import server as gmsh_server

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
        (gmsh_server, "gmsh_audit_summary"),
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
        assert summary["top_files"][0]["path"] == "examples\\a.py"


def test_gmsh_numsubedges_remediation_plan(monkeypatch, tmp_path):
    from radia_mcp.gmsh import server

    examples = tmp_path / "examples"
    examples.mkdir()
    target = examples / "curved.py"
    target.write_text("mesh.Curve(3)\n", encoding="utf-8")
    clean = examples / "flat.py"
    clean.write_text("print('flat')\n", encoding="utf-8")

    def fake_lint(filepath: str):
        if filepath.endswith("curved.py"):
            return [{"line": 0, "severity": "MODERATE",
                     "rule": "numsubedges-missing", "message": "x"}]
        return []

    monkeypatch.setattr(server, "PROJECT_ROOT", tmp_path)
    monkeypatch.setattr(server, "_lint_file", fake_lint)

    plan = server.gmsh_numsubedges_remediation_plan("examples", limit=1)
    assert plan["ok"] is True
    assert plan["total_affected"] == 1
    assert plan["returned"] == 1
    assert plan["truncated"] is False
    assert plan["directory_groups"] == [{
        "directory": "examples",
        "count": 1,
        "directory_companion": {
            "geo_companion": "examples\\_gmsh_display.geo",
            "geo_template": (
                "// Shared GMSH display companion for examples\n"
                "// Use with any high-order .msh output from this directory.\n"
                "Mesh.NumSubEdges = 4;\n"
                "// Merge \"<result>.msh\";\n"
            ),
        },
    }]
    item = plan["affected"][0]
    assert item["script"] == "examples\\curved.py"
    assert item["triggers"] == ["high_order_curve"]
    assert item["geo_companion"] == "examples\\curved_display.geo"
    assert "Mesh.NumSubEdges = 4;" in item["geo_template"]


def test_gmsh_numsubedges_rule_respects_display_companion(tmp_path):
    from radia_mcp.gmsh.rules import check_numsubedges_missing

    script = tmp_path / "curved.py"
    lines = ["mesh.Curve(3)\n"]
    assert check_numsubedges_missing(str(script), lines)

    (tmp_path / "curved_display.geo").write_text(
        "Mesh.NumSubEdges = 4;\n",
        encoding="utf-8",
    )
    assert check_numsubedges_missing(str(script), lines) == []


def test_gmsh_numsubedges_rule_respects_directory_display_companion(tmp_path):
    from radia_mcp.gmsh.rules import check_numsubedges_missing

    script = tmp_path / "curved.py"
    lines = ["mesh.Curve(3)\n"]
    assert check_numsubedges_missing(str(script), lines)

    (tmp_path / "_gmsh_display.geo").write_text(
        "Mesh.NumSubEdges = 4;\n",
        encoding="utf-8",
    )
    assert check_numsubedges_missing(str(script), lines) == []


def test_gmsh_numsubedges_rule_respects_ancestor_display_companion(tmp_path):
    from radia_mcp.gmsh.rules import check_numsubedges_missing

    examples = tmp_path / "examples"
    nested = examples / "nested" / "case"
    nested.mkdir(parents=True)
    script = nested / "curved.py"
    lines = ["mesh.Curve(3)\n"]
    assert check_numsubedges_missing(str(script), lines)

    (examples / "_gmsh_display.geo").write_text(
        "Mesh.NumSubEdges = 4;\n",
        encoding="utf-8",
    )
    assert check_numsubedges_missing(str(script), lines) == []


def test_gmsh_mesh_generation_remediation_plan(monkeypatch, tmp_path):
    from radia_mcp.gmsh import server

    examples = tmp_path / "examples"
    examples.mkdir()
    target = examples / "makes_mesh.py"
    target.write_text(
        "import gmsh\n"
        "gmsh.model.occ.addBox(0, 0, 0, 1, 1, 1)\n"
        "gmsh.model.mesh.generate(3)\n",
        encoding="utf-8",
    )
    clean = examples / "display_only.py"
    clean.write_text("print('display')\n", encoding="utf-8")

    def fake_lint(filepath: str):
        if filepath.endswith("makes_mesh.py"):
            return [
                {"line": 2, "severity": "CRITICAL",
                 "rule": "gmsh-mesh-generation", "message": "occ"},
                {"line": 3, "severity": "CRITICAL",
                 "rule": "gmsh-mesh-generation", "message": "mesh"},
            ]
        return []

    monkeypatch.setattr(server, "PROJECT_ROOT", tmp_path)
    monkeypatch.setattr(server, "_lint_file", fake_lint)

    plan = server.gmsh_mesh_generation_remediation_plan("examples", limit=1)
    assert plan["ok"] is True
    assert plan["total_affected"] == 1
    assert plan["total_findings"] == 2
    assert plan["directory_groups"] == [{"directory": "examples", "findings": 2}]
    item = plan["affected"][0]
    assert item["script"] == "examples\\makes_mesh.py"
    assert item["findings"][0]["line"] == 2
    assert "gmsh.model.occ.addBox" in item["findings"][0]["snippet"]
    assert "Mesh('model.vol')" in item["mesh_output_hint"]
