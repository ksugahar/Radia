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
        assert set(triggers[event]["paths"]) == set(IMPACT.NATIVE_PATHS) | {
            "pyproject.toml", "tests/test_tier_manifest.json"}
    steps = workflow["jobs"]["ams-regression"]["steps"]
    selector = next(s for s in steps if s.get("id") == "native-impact")
    assert "python tools/sparsesolv_ci_impact.py" in selector["run"]
    for step in steps:
        if step.get("name", "").startswith(("Create isolated", "Build ", "Verify source", "Upload solver")):
            assert "steps.native-impact.outputs.required == 'true'" in step["if"]


def manifest():
    return {"schema": "radia.test-tier-manifest.v1", "impact_rules": {}, "profiles": {
        "sparsesolv": {"extends": "base", "paths": [SUITE], "max_elapsed_seconds": 30},
        "base": {"paths": ["tests/base.py"]}, "fast": {"paths": ["tests/old.py"]}}}


@pytest.mark.parametrize("change,expected", [
    (lambda d: d["impact_rules"].update({"tools/new.py": ["tests/new.py"]}), True),
    (lambda d: d["profiles"]["fast"]["paths"].append("tests/new.py"), True),
    (lambda d: d["profiles"]["sparsesolv"]["paths"].append("tests/new.py"), False),
    (lambda d: d["profiles"]["base"]["paths"].append("tests/new.py"), False),
    (lambda d: d["profiles"]["sparsesolv"].update(max_elapsed_seconds=60), False),
    (lambda d: d["profiles"]["sparsesolv"].update(extends="fast"), False),
    (lambda d: d.update(future_global_option=True), False),
])
def test_manifest_compares_ams_and_all_ancestors(change, expected):
    before, after = manifest(), manifest()
    change(after)
    assert IMPACT.ams_manifest_unchanged(json.dumps(before), json.dumps(after)) is expected


@pytest.mark.parametrize("change", [
    lambda d: d.update(schema="unknown"),
    lambda d: d["profiles"].pop("sparsesolv"),
    lambda d: d["profiles"].pop("base"),
    lambda d: d["profiles"]["base"].update(extends="sparsesolv"),
    lambda d: d["profiles"]["base"].update(paths="bad"),
    lambda d: d["profiles"]["base"].update(paths=[SUITE]),
    lambda d: d["profiles"]["sparsesolv"].update(max_elapsed_seconds=-1),
])
def test_invalid_manifest_fails_closed(change):
    before, after = manifest(), manifest()
    change(after)
    with pytest.raises((KeyError, ValueError, TypeError)):
        IMPACT.ams_manifest_unchanged(json.dumps(before), json.dumps(after))


@pytest.mark.parametrize("extra,required", [
    ("tests/new.py\0tools/new.py\0", False),
    ("src/ext/sparsesolv/test.cpp\0", True),
    ("matlab/+radia/+internal/callMex.m\0", True),
    ("tools/run_test_tier.py\0", True),
    (".github/workflows/sparsesolv.yml\0", True),
])
def test_manifest_edit_with_other_files(extra, required):
    path = "tests/test_tier_manifest.json"
    before, after = manifest(), manifest()
    after["impact_rules"]["tools/new.py"] = ["tests/new.py"]
    def git(*args):
        if args[0] == "diff":
            return path + "\0" + extra
        return json.dumps(before if args[1].startswith("base:") else after)
    assert IMPACT.native_required("push", {"before": "base"}, git)[0] is required


def test_mixed_manifest_and_runtime_dependency_change_still_runs():
    def git(*args):
        if args[0] == "diff":
            return "tests/test_tier_manifest.json\0pyproject.toml\0tests/new.py\0"
        if args[1].endswith("pyproject.toml"):
            return PROJECT if args[1].startswith("base:") else PROJECT.replace("6.2.2606", "6.2.2607")
        return json.dumps(manifest())
    assert IMPACT.native_required("push", {"before": "base"}, git)[0]
