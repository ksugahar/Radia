# -*- coding: utf-8 -*-
"""Path-guard tests for radia-ngsolve MCP lint tools."""

import os
import sys


_SRC = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "src"))
if _SRC not in sys.path:
    sys.path.insert(0, _SRC)

from radia_mcp.radia_ngsolve.server import (  # noqa: E402
    lint_radia_directory,
    lint_radia_script,
)


def test_lint_radia_script_rejects_unregistered_absolute_path(tmp_path):
    secret = tmp_path / "secret.py"
    secret.write_text("print('not a public lint root')\n", encoding="utf-8")

    result = lint_radia_script(str(secret))

    assert result.startswith("Error: Access denied:")
    assert "RADIA_MCP_LINT_ROOTS" in result
    assert str(tmp_path) not in result
    assert "<outside allowed lint roots>" in result


def test_lint_radia_script_allows_explicit_opt_in_root(tmp_path, monkeypatch):
    script = tmp_path / "case.py"
    script.write_text("print('ok')\n", encoding="utf-8")
    monkeypatch.setenv("RADIA_MCP_LINT_ROOTS", str(tmp_path))

    result = lint_radia_script(str(script))

    assert result.startswith("[OK]")
    assert "case.py" in result
    assert str(tmp_path) not in result


def test_lint_radia_directory_rejects_unregistered_absolute_path(tmp_path):
    (tmp_path / "case.py").write_text("print('not allowed')\n", encoding="utf-8")

    result = lint_radia_directory(str(tmp_path))

    assert result.startswith("Error: Access denied:")


def test_lint_radia_script_accepts_monorepo_relative_package_path():
    result = lint_radia_script(
        "packages/radia-mcp/src/radia_mcp/radia_ngsolve/server.py"
    )

    assert not result.startswith("Error: File not found:")
    assert "repo:/packages/radia-mcp/src/radia_mcp/radia_ngsolve/server.py" in result
    assert ":\\" not in result.splitlines()[0]
