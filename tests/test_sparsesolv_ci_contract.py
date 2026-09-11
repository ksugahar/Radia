"""Keep native solver coverage connected without importing native dependencies."""

import json
import importlib.util
from pathlib import Path

import pytest
import yaml

ROOT = Path(__file__).resolve().parents[1]
SUITE = "src/ext/sparsesolv/tests/test_sparsesolv.py"
SPEC = importlib.util.spec_from_file_location("ams_impact", ROOT / "tools/sparsesolv_ci_impact.py")
IMPACT = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(IMPACT)
PROJECT = '''[build-system]
requires = ["setuptools"]
[project]
dependencies = ["ngsolve==6.2.2606"]
[project.optional-dependencies]
test = ["pytest"]
dev = ["ruff"]
beam = ["scipy"]
'''


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


@pytest.mark.parametrize("old,new,skip", [
    ('test = ["pytest"]', 'test = ["pytest", "setuptools"]', True),
    ('dev = ["ruff"]', 'dev = []', True),
    ('beam = ["scipy"]', 'beam = []', False),
    ('6.2.2606', '6.2.2607', False),
    ('requires = ["setuptools"]', 'requires = ["setuptools>=80"]', False),
    ('test = ["pytest"]', 'test = "invalid"', False),
    ('test = ["pytest"]', 'test = ["pytest"] # comment', False),
])
def test_native_impact_only_exempts_test_and_dev_extras(old, new, skip):
    assert IMPACT.test_dependencies_only(PROJECT, PROJECT.replace(old, new)) is skip


def test_native_impact_preserves_mixed_configuration_changes():
    after = PROJECT.replace('["pytest"]', '[]').replace('6.2.2606', '6.2.2607')
    assert not IMPACT.test_dependencies_only(PROJECT, after)


@pytest.mark.parametrize("event,base", [("push", "a" * 40), ("pull_request", "HEAD^1")])
def test_native_impact_uses_complete_event_range(event, base):
    calls = []
    def git(*args):
        calls.append(args)
        return {("diff", "--name-only", "-z", base, "HEAD", "--"): "pyproject.toml\0",
                ("show", f"{base}:pyproject.toml"): PROJECT,
                ("show", "HEAD:pyproject.toml"): PROJECT.replace('["pytest"]', '[]')}[args]
    assert IMPACT.native_required(event, {"before": base}, git)[0] is False
    assert calls[0][3] == base


@pytest.mark.parametrize("paths", ["", "pyproject.toml\0Build.ps1\0", "src/ext/sparsesolv/test.cpp\0"])
def test_native_impact_does_not_exempt_other_files(paths):
    assert IMPACT.native_required("push", {"before": "a" * 40}, lambda *args: paths)[0]


def test_manual_native_impact_always_runs():
    assert IMPACT.native_required("workflow_dispatch", {}, None)[0]
    assert IMPACT.native_required("push", {"before": "0" * 40}, None)[0]


@pytest.mark.parametrize("error", [ValueError("invalid TOML"), OSError("missing git base")])
def test_unknown_native_impact_fails_closed(tmp_path, monkeypatch, error):
    event = tmp_path / "event.json"
    event.write_text('{"before":"base"}')
    output = tmp_path / "output"
    monkeypatch.setenv("GITHUB_EVENT_PATH", str(event))
    monkeypatch.setenv("GITHUB_EVENT_NAME", "push")
    monkeypatch.setenv("GITHUB_OUTPUT", str(output))
    def fail(*args):
        raise error
    monkeypatch.setattr(IMPACT, "native_required", fail)
    IMPACT.main()
    assert output.read_text() == "required=true\n"


def test_native_impact_gates_expensive_steps_and_watches_itself():
    workflow = yaml.safe_load((ROOT / ".github/workflows/sparsesolv.yml").read_text())
    triggers = workflow.get("on", workflow.get(True))
    for event in ("push", "pull_request"):
        assert "tools/sparsesolv_ci_impact.py" in triggers[event]["paths"]
    steps = workflow["jobs"]["ams-regression"]["steps"]
    selector = next(s for s in steps if s.get("id") == "native-impact")
    assert "python tools/sparsesolv_ci_impact.py" in selector["run"]
    for step in steps:
        if step.get("name", "").startswith(("Create isolated", "Build ", "Verify source", "Upload solver")):
            assert "steps.native-impact.outputs.required == 'true'" in step["if"]
