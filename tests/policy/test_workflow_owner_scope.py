"""A new workflow must explicitly declare whether it belongs to solver CI."""
from fnmatch import fnmatchcase
from pathlib import Path

import yaml

ROOT = Path(__file__).resolve().parents[2]
SOLVER = {"radia-fast.yml", "build-test.yml", "release.yml"}
INDEPENDENT = {
    "eqnedit64.yml", "release-eqnedit64-pypi.yml", "radia-mcp-matrix.yml",
    "cubit-mesh-export.yml", "radia-optuna.yml", "release-radia-optuna.yml",
    "radia-gpu.yml", "radia-gpu-release.yml", "sparsesolv.yml",
    "policy-lint.yml", "ltspice-data-safety.yml",
}


def test_every_workflow_has_an_explicit_owner():
    names = {p.name for p in (ROOT / ".github/workflows").glob("*.yml")}
    assert names == SOLVER | INDEPENDENT
    assert not SOLVER & INDEPENDENT


def test_solver_triggers_exclude_only_independent_workflows_and_policy_tests():
    workflow = yaml.load((ROOT / ".github/workflows/radia-fast.yml").read_text(), Loader=yaml.BaseLoader)
    for event in ("push", "pull_request"):
        ignored = workflow["on"][event]["paths-ignore"]
        excluded = lambda path: any(fnmatchcase(path, pattern) for pattern in ignored)
        assert all(excluded(".github/workflows/" + name) for name in INDEPENDENT)
        assert all(not excluded(".github/workflows/" + name) for name in SOLVER)
        assert excluded("tests/policy/test_release_workflow_ref_gate.py")
        assert excluded("packages/eqnedit64/src/eqnedit64/api.py")
        assert excluded("tools/eqnedit64/src/office_html.cpp")
        assert excluded("packages/radia-mcp/src/radia_mcp/presentation/tools.py")
        for path in ("src/radia/__init__.py", "Build.ps1", "tests/test_wheel_package_policy.py", "tools/ci_preflight.py"):
            assert not excluded(path), path


def test_policy_contracts_run_in_policy_lane():
    workflow = (ROOT / ".github/workflows/policy-lint.yml").read_text()
    assert "python -m pytest tests/policy/ -q" in workflow
    manifest = (ROOT / "tests/test_tier_manifest.json").read_text()
    assert "test_release_workflow_ref_gate.py" not in manifest
