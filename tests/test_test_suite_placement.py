"""Structural guard for the short-tests/heavy-validation boundary."""

from __future__ import annotations

import ast
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def _slow_locations(path: Path) -> list[int]:
    tree = ast.parse(path.read_text(encoding="utf-8"), filename=str(path))
    locations: list[int] = []
    for node in ast.walk(tree):
        if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef, ast.ClassDef)):
            if not node.name.startswith("test_") and not node.name.startswith("Test"):
                continue
            if any(ast.unparse(item).endswith(".slow") for item in node.decorator_list):
                locations.append(node.lineno)
        elif isinstance(node, ast.Assign):
            if any(isinstance(target, ast.Name) and target.id == "pytestmark"
                   for target in node.targets):
                if ".slow" in ast.unparse(node.value):
                    locations.append(node.lineno)
    return locations


def test_short_suite_does_not_hide_slow_tests():
    violations = {}
    for path in (ROOT / "tests").rglob("test_*.py"):
        locations = _slow_locations(path)
        if locations:
            violations[str(path.relative_to(ROOT))] = locations
    assert violations == {}
    assert not (ROOT / "tests" / "ci_slow_nodeids.txt").exists()


def test_measured_slow_nodeids_belong_to_validation_suite():
    path = ROOT / "validation_test" / "slow_nodeids.txt"
    nodeids = [
        line.split("#", 1)[0].strip().replace("\\", "/")
        for line in path.read_text(encoding="utf-8").splitlines()
        if line.split("#", 1)[0].strip()
    ]
    assert nodeids
    assert all(nodeid.startswith("validation_test/") for nodeid in nodeids)
