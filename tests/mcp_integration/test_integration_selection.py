"""Changed consumed inputs must stay in the change-scoped CI plan.

Guards: Missing integration checks after foreign-input moves, and unrelated
foreign changes expanding the MCP lane.
Assumes: Source selection only; workflow/owner projections have separate audits.
"""

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
    "tests/mcp_integration/test_matlab_source_contract.py",
    "validation_test/radia_mcp/artifacts/annular_motor_dual_lane_v1/manifest.json",
    "validation_test/vim_coupled/results_disk.json",
    "packages/radia-simulink/matlab/+radia/+acoustic/cqGrid.m",
    "packages/radia-optuna/matlab/optuna50_api_coverage.json",
    "packages/radia-optuna/oracles/optuna50_oracle.json",
    "packages/radia-optuna/oracles/generate_optuna50_mcp_oracle.py",
    "packages/radia-optuna/tools/generate_optuna50_api_coverage.py",
    "tests/matlab/test_optuna_sampler.m",
    "tests/mcp_server/fixtures/bad_radia_script.py",
    "validation_test/radia_mcp/validation_vol_tet_quality.py",
    "packages/radia-mcp/src/radia_mcp/radia_ngsolve/netgen_vol.py",
    "packages/radia-mcp/src/radia_mcp/radia_ngsolve/force.py",
    "packages/radia-mcp/src/radia_mcp/radia_ngsolve/scalar_fem3d.py",
])
def test_repository_changes_keep_integration_lane(changed):
    assert SELECTOR.build_plan([changed])["integration_tests"] == ["tests/mcp_integration"]


@pytest.mark.parametrize("changed", [
    "validation_test/maglev/results.json",
    "packages/radia-simulink/matlab/+radia/contract.m",
    "packages/radia-optuna/src/radia_optuna/__init__.py",
    "packages/radia-optuna/tools/build_mexcuda.m",
    "packages/radia-optuna/oracles/README.md",
    "docs/README.md",
])
def test_unread_foreign_changes_do_not_run_the_integration_lane(changed):
    assert SELECTOR.build_plan([changed])["integration_tests"] == []


def test_domain_changes_keep_moved_contracts():
    plan = SELECTOR.build_plan(["packages/radia-mcp/src/radia_mcp/matlab/server.py"])
    assert "tests/mcp_integration/test_matlab_source_contract.py" in plan["integration_tests"]
    assert "tests/mcp_integration/test_matlab_optuna_quality.py" in plan["integration_tests"]
    assert "tests/mcp_integration/test_cubit_golden.py" not in plan["integration_tests"]


def test_full_audit_includes_integration():
    assert SELECTOR.build_plan([], full=True)["integration_tests"] == ["tests/mcp_integration"]
