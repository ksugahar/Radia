"""Radia-only test lanes must not import exporter-owned Cubit modules."""

import ast
from pathlib import Path


def test_radia_test_lanes_do_not_import_cubit_exporter():
    root = Path(__file__).resolve().parents[2]
    violations = []
    for directory in ("packages/radia-mcp/tests", "tests/mcp_server", "tests/mcp_integration"):
        for path in (root / directory).rglob("test_*.py"):
            for node in ast.walk(ast.parse(path.read_text(encoding="utf-8-sig"))):
                modules = ([node.module or ""] if isinstance(node, ast.ImportFrom)
                           else [a.name for a in node.names] if isinstance(node, ast.Import)
                           else [])
                if any(m == "cubit_mesh_export" or m.startswith("cubit_mesh_export.")
                       for m in modules):
                    violations.append(f"{path.relative_to(root)}:{node.lineno}")
    assert not violations, "Move Cubit-owned contracts to exporter tests: " + str(violations)
