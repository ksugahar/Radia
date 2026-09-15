"""Moved tests must stay in the change-scoped CI plan."""

import importlib.util
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[2]
SPEC = importlib.util.spec_from_file_location(
    "integration_selector", ROOT / "packages/radia-mcp/tools/select_ci_tests.py"
)
SELECTOR = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(SELECTOR)


@pytest.mark.parametrize("changed", [
    "tests/mcp_integration/test_cubit_golden.py",
    "validation_test/maglev/results.json",
    "matlab/+radia/contract.m",
    "tests/mcp_server/fixtures/bad_cubit_script.py",
])
def test_repository_changes_keep_integration_lane(changed):
    assert SELECTOR.build_plan([changed])["integration_tests"] == ["tests/mcp_integration"]


def test_domain_changes_keep_moved_contracts():
    plan = SELECTOR.build_plan(["packages/radia-mcp/src/radia_mcp/matlab/server.py"])
    assert "tests/mcp_integration/test_matlab_source_contract.py" in plan["integration_tests"]
    assert "tests/mcp_integration/test_matlab_optuna_quality.py" in plan["integration_tests"]
    assert "tests/mcp_integration/test_cubit_golden.py" not in plan["integration_tests"]


def test_full_audit_includes_integration():
    assert SELECTOR.build_plan([], full=True)["integration_tests"] == ["tests/mcp_integration"]
