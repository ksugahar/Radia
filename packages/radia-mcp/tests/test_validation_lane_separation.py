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
def test_package_tests_do_not_walk_above_the_package():
    """Direct monorepo traversal belongs in tests/mcp_integration instead."""
    violations = []
    for path in Path(__file__).parent.rglob("test_*.py"):
        tree = ast.parse(path.read_text(encoding="utf-8"))
        for node in ast.walk(tree):
            if not isinstance(node, ast.Subscript):
                continue
            value = node.value
            if (isinstance(value, ast.Attribute) and value.attr == "parents"
                    and "__file__" in ast.unparse(value)
                    and isinstance(node.slice, ast.Constant)
                    and isinstance(node.slice.value, int) and node.slice.value >= 2):
                violations.append(f"{path.name}:{node.lineno}")
    assert not violations, f"Package tests escape their source tree: {violations}"
