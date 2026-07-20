from __future__ import annotations

import ast
from pathlib import Path


def _solver_imports(path: Path) -> set[str]:
    tree = ast.parse(path.read_text(encoding="utf-8", errors="ignore"))
    imports: set[str] = set()
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            imports.update(alias.name.split(".", 1)[0] for alias in node.names)
        elif isinstance(node, ast.ImportFrom) and node.module:
            imports.add(node.module.split(".", 1)[0])
        elif (
            isinstance(node, ast.Call)
            and isinstance(node.func, ast.Attribute)
            and node.func.attr == "importorskip"
            and node.args
            and isinstance(node.args[0], ast.Constant)
        ):
            imports.add(str(node.args[0].value).split(".", 1)[0])
    return imports & {"netgen", "ngsolve"}


def test_solver_backed_numerical_checks_use_validation_lane():
    test_root = Path(__file__).resolve().parent
    offenders = {}
    for path in test_root.rglob("test_*.py"):
        imports = _solver_imports(path)
        if imports:
            offenders[path.name] = sorted(imports)

    assert offenders == {}, (
        "Move solver-backed numerical checks to validation_test/radia_mcp; "
        f"package tests must remain fast API/MCP contracts: {offenders}"
    )
