"""Contracts for impact-scoped radia-mcp CI selection."""

from __future__ import annotations

import importlib.util
import json
from pathlib import Path
import subprocess
import sys


PACKAGE_ROOT = Path(__file__).resolve().parents[1]
SELECTOR_PATH = PACKAGE_ROOT / "tools" / "select_ci_tests.py"
WORKFLOW_PATH = (
    PACKAGE_ROOT.parents[1] / ".github" / "workflows" / "radia-mcp-matrix.yml"
)
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
    assert "tests/test_common_mcp_runtime_contract.py" in plan["package_tests"]
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


def test_server_registration_change_does_not_select_numerical_family_suite():
    plan = SELECTOR.build_plan(
        ["packages/radia-mcp/src/radia_mcp/build123d/server.py"]
    )
    assert plan["server_selftests"] == ["build123d"]
    assert "tests/test_build123d_operations.py" not in plan["package_tests"]
    assert len(plan["package_tests"]) < 60


def test_large_server_registration_change_stays_on_contract_lane():
    plan = SELECTOR.build_plan(
        ["packages/radia-mcp/src/radia_mcp/radia_ngsolve/server.py"]
    )
    assert plan["server_selftests"] == ["radia-ngsolve"]
    assert "tests/test_radia_ngsolve_mcp_contract.py" in plan["package_tests"]
    assert "tests/test_force_coenergy_gate.py" not in plan["package_tests"]
    assert len(plan["package_tests"]) < 20


def test_common_change_selftests_every_server_without_selecting_every_test():
    plan = SELECTOR.build_plan(
        ["packages/radia-mcp/src/radia_mcp/common/status.py"]
    )
    assert len(plan["server_selftests"]) >= 30
    assert "meta" in plan["server_selftests"]
    assert any(
        selector.startswith("tests/test_meta_health.py::")
        for selector in plan["package_tests"]
    )
    assert len(plan["package_tests"]) < len(list((PACKAGE_ROOT / "tests").glob("test_*.py")))
    assert not any(
        "::" in selector
        and selector.split("::", 1)[0] in plan["package_tests"]
        for selector in plan["package_tests"]
    )


def test_common_lazy_export_change_skips_unrelated_rag_content_tests():
    plan = SELECTOR.build_plan(
        ["packages/radia-mcp/src/radia_mcp/common/__init__.py"]
    )
    assert "tests/test_coarse_tool_registry.py" in plan["package_tests"]
    assert "tests/test_optional_dependency_imports.py" in plan["package_tests"]
    assert "tests/test_chroma_multilingual.py" not in plan["package_tests"]
    assert "tests/test_chunk_garble_gate.py" not in plan["package_tests"]


def test_changed_symbol_narrows_a_large_compatibility_module():
    path = "packages/radia-mcp/src/radia_mcp/radia_ngsolve/solve.py"
    plan = SELECTOR.build_plan(
        [path],
        changed_symbols_by_file={path: {"laminated_mu_eff"}},
    )
    assert "tests/test_lamination_adapter.py" in plan["package_tests"]
    assert "tests/test_two_wire_force.py" not in plan["package_tests"]
    assert len(plan["package_tests"]) < 15


def test_changed_lines_map_to_function_or_module_scope():
    source = "value = 1\n\ndef alpha():\n    return 1\n\ndef beta():\n    return 2\n"
    assert SELECTOR._symbols_for_line_ranges(source, [(4, 4)]) == {"alpha"}
    assert SELECTOR._symbols_for_line_ranges(source, [(1, 1)]) == {"*"}


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


def test_cli_accepts_large_changed_file_list_over_stdin():
    changed = [
        f"packages/radia-mcp/src/radia_mcp/family_{index}/server.py"
        for index in range(2000)
    ]
    completed = subprocess.run(
        [sys.executable, str(SELECTOR_PATH), "--changed-files-json", "-"],
        input=json.dumps(changed),
        capture_output=True,
        text=True,
        check=True,
    )
    plan = json.loads(completed.stdout)
    assert plan["mode"] == "targeted"
    assert plan["changed_files"] == sorted(changed)


def test_workflow_discovers_selected_tests_from_the_test_directory():
    workflow = WORKFLOW_PATH.read_text(encoding="utf-8")

    assert 'env["RADIA_MCP_CI_SELECTION_JSON"]' in workflow
    assert '"-m", "not xval and not slow", "tests"' in workflow
    assert '*targets' not in workflow
