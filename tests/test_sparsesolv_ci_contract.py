"""Keep native solver coverage connected without importing native dependencies."""

import json
from pathlib import Path

import yaml

ROOT = Path(__file__).resolve().parents[1]
SUITE = "src/ext/sparsesolv/tests/test_sparsesolv.py"


def test_sparsesolv_is_in_native_release_and_focused_tiers_only():
    profiles = json.loads((ROOT / "tests/test_tier_manifest.json").read_text())["profiles"]
    assert SUITE in profiles["native-smoke"]["paths"]
    assert profiles["sparsesolv"]["paths"] == [SUITE]
    assert SUITE not in profiles["fast-contracts"]["paths"]


def test_solver_workflow_builds_before_testing_on_mdx():
    workflow = yaml.safe_load((ROOT / ".github/workflows/sparsesolv.yml").read_text())
    triggers = workflow.get("on", workflow.get(True))
    for event in ("push", "pull_request"):
        assert "src/ext/sparsesolv/**" in triggers[event]["paths"]
        assert "pyproject.toml" in triggers[event]["paths"]
    job = workflow["jobs"]["ams-regression"]
    assert "mdx" in job["runs-on"]
    assert "head.repo.full_name == github.repository" in job["if"]
    script = "\n".join(step.get("run", "") for step in job["steps"])
    assert script.index("--target sparsesolv_ngsolve") < script.index("--profile sparsesolv")
    assert "s.has_compact_ams()" in script
    assert "p.parent ==" in script


def test_ams_unavailable_is_not_a_silent_skip():
    source = (ROOT / SUITE).read_text()
    assert "pytest.skip(" not in source
    assert source.count("assert has_compact_ams()") == 3
