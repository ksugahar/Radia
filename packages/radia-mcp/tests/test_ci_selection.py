"""Contracts for impact-scoped radia-mcp CI selection."""

from __future__ import annotations

import importlib.util
from pathlib import Path


PACKAGE_ROOT = Path(__file__).resolve().parents[1]
SELECTOR_PATH = PACKAGE_ROOT / "tools" / "select_ci_tests.py"
SPEC = importlib.util.spec_from_file_location("select_ci_tests", SELECTOR_PATH)
assert SPEC and SPEC.loader
SELECTOR = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(SELECTOR)


def test_contract_set_is_always_selected():
    plan = SELECTOR.build_plan(["packages/radia-mcp/README.md"])
    assert plan["mode"] == "targeted"
    assert "tests/test_ci_selection.py" in plan["package_tests"]
    assert any(
        selector.startswith("tests/test_meta_health.py::")
        for selector in plan["package_tests"]
    )
    assert "tests/test_meta_health.py" not in plan["package_tests"]
    assert plan["server_selftests"] == []


def test_changed_test_file_is_selected_directly():
    path = "packages/radia-mcp/tests/test_force_mcp_contract.py"
    plan = SELECTOR.build_plan([path])
    assert "tests/test_force_mcp_contract.py" in plan["package_tests"]
    assert "tests/test_validation_lane_separation.py" in plan["package_tests"]


def test_server_family_change_selects_related_tests_and_one_server():
    plan = SELECTOR.build_plan(
        ["packages/radia-mcp/src/radia_mcp/build123d/operations.py"]
    )
    assert plan["server_selftests"] == ["build123d"]
    assert any("build123d" in path for path in plan["package_tests"])
    assert "tests" not in plan["package_tests"]


def test_common_change_selftests_every_server_without_selecting_every_test():
    plan = SELECTOR.build_plan(
        ["packages/radia-mcp/src/radia_mcp/common/status.py"]
    )
    assert len(plan["server_selftests"]) >= 30
    assert "meta" in plan["server_selftests"]
    assert "tests/test_meta_health.py" in plan["package_tests"]
    assert len(plan["package_tests"]) < len(list((PACKAGE_ROOT / "tests").glob("test_*.py")))
    assert not any(
        "::" in selector
        and selector.split("::", 1)[0] in plan["package_tests"]
        for selector in plan["package_tests"]
    )


def test_package_metadata_selects_packaging_contracts_only():
    plan = SELECTOR.build_plan(["packages/radia-mcp/pyproject.toml"])
    assert "tests/test_mcp_sdk_dependency_contract.py" in plan["package_tests"]
    assert "tests/test_optional_dependency_imports.py" in plan["package_tests"]
    assert plan["server_selftests"] == []


def test_explicit_full_audit_selects_whole_suite_and_catalog():
    plan = SELECTOR.build_plan([], full=True)
    assert plan["mode"] == "full"
    assert plan["package_tests"] == ["tests"]
    assert len(plan["server_selftests"]) >= 30
    assert plan["run_mcp_response_tests"] is True


def test_root_mcp_response_change_enables_only_its_regression_lane():
    plan = SELECTOR.build_plan(["tests/mcp_server/test_mcp_tool_responses.py"])
    assert plan["run_mcp_response_tests"] is True
    assert plan["server_selftests"] == []
