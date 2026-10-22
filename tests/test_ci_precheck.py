"""Precheck coverage must match CI and propagate real selected failures."""
import importlib.util
import json
from pathlib import Path
import shlex
import subprocess
import sys

import pytest
import yaml

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location("source_precheck", ROOT / "tools/ci_precheck.py")
precheck = importlib.util.module_from_spec(spec)
spec.loader.exec_module(precheck)


def git_diff(paths):
    def git(*args):
        if args[0] == "merge-base":
            return "a" * 40
        assert "--no-renames" in args and "-z" in args
        return "\0".join(paths) + ("\0" if paths else "")
    return git


def test_policy_collection_is_the_workflows_actual_collection():
    workflow = yaml.load((ROOT / ".github/workflows/policy-lint.yml").read_text(), Loader=yaml.BaseLoader)
    commands = [line.strip() for step in workflow["jobs"]["policy-checks"]["steps"]
                for line in step.get("run", "").splitlines()]
    pytest_commands = [shlex.split(line) for line in commands if line.startswith("python -m pytest ")]
    assert len(pytest_commands) == 1
    targets = tuple(word for word in pytest_commands[0][3:] if not word.startswith("-"))
    assert tuple(path.rstrip("/") for path in precheck.POLICY_TESTS) == tuple(path.rstrip("/") for path in targets)
    assert "python tools/policy_lint.py" in commands
    assert "python tools/eqnedit64/tests/check_public_names.py" in commands


def test_docs_dont_create_an_unneeded_light_environment():
    plan = precheck.build_source_plan("base", git=git_diff(["packages/cubit-mesh-export/README.md"]))
    assert not plan["radia"] and not plan["mcp"]


@pytest.mark.parametrize("path,lane", [
    ("src/radia/bh_law.py", "radia"),
    ("packages/radia-mcp/src/radia_mcp/bibliography/metadata.py", "mcp"),
    ("tests/mcp_integration/test_netgen_vol_frozen_regression.py", "mcp"),
])
def test_real_distribution_inputs_select_their_source_checks(path, lane):
    assert precheck.build_source_plan("base", git=git_diff([path]))[lane]


def test_shared_precheck_changes_cannot_skip_the_precheck_or_mcp():
    plan = precheck.build_source_plan("base", git=git_diff(["tools/ci_precheck.py"]))
    assert plan["radia"] and plan["mcp"]


def test_rename_delete_diff_retains_both_sides():
    plan = precheck.build_source_plan("base", git=git_diff(["src/radia/bh_law.py", "docs/moved.md"]))
    assert plan["changed_paths"] == ["docs/moved.md", "src/radia/bh_law.py"]
    assert plan["radia"]


def test_unknown_ownership_cannot_become_a_successful_skip():
    with pytest.raises(ValueError):
        precheck.build_source_plan("base", git=git_diff(["unknown_new_owner/input.dat"]))


def test_full_does_not_require_a_git_comparison_and_empty_diff_keeps_source_baseline():
    def unavailable(*args):
        raise AssertionError("full must not need a base")
    assert precheck.build_source_plan("missing", full=True, git=unavailable)["mcp"]
    assert precheck.build_source_plan("base", git=git_diff([]))["radia"]


def fake_run(commands, *, selection=None, failing=None):
    def run(name, command, env, steps, cwd=None):
        commands.append((name, command, env, cwd))
        if name == "MCP test selection":
            Path(command[command.index("--output") + 1]).write_text(json.dumps(selection))
        rc = int(name == failing)
        steps.append({"name": name, "returncode": rc})
        return not rc
    return run


def test_policy_always_runs_and_plan_only_never_claims_tests_passed(tmp_path, monkeypatch):
    plan = precheck.build_source_plan("base", git=git_diff(["packages/cubit-mesh-export/README.md"]))
    monkeypatch.setattr(precheck, "build_source_plan", lambda *a, **k: plan)
    commands = []
    monkeypatch.setattr(precheck, "_run", fake_run(commands))
    report = tmp_path / "report.json"
    assert precheck.main(["--report", str(report)]) == 0
    assert [row[0] for row in commands] == ["policy lint", "public-name lint", "policy tests", "repository contracts"]
    assert precheck.main(["--plan-only", "--report", str(report)]) == 0
    result = json.loads(report.read_text())
    assert result["status"] == "planned" and result["passed"] is None and not result["qualified"]


def test_invalid_base_records_failure_before_any_execution(tmp_path, monkeypatch):
    def unavailable(*args, **kwargs):
        raise subprocess.CalledProcessError(128, ["git", "merge-base"])
    monkeypatch.setattr(precheck, "build_source_plan", unavailable)
    monkeypatch.setattr(precheck, "_run", lambda *a, **k: pytest.fail("must not execute"))
    report = tmp_path / "report.json"
    assert precheck.main(["--report", str(report)]) == 1
    assert not json.loads(report.read_text())["passed"]


def selection(**changes):
    result = {"schema": "radia-mcp.ci-selection.v1",
              "package_tests": ["tests/test_ci_selection.py"],
              "integration_tests": ["tests/mcp_integration/test_contract.py"],
              "server_selftests": [], "run_mcp_response_tests": True}
    result.update(changes)
    return result


def test_mcp_executes_the_actual_collection_protocol_and_propagates_boundary_failure(tmp_path, monkeypatch):
    commands = []
    monkeypatch.setattr(precheck, "_run", fake_run(commands, selection=selection(),
                                                 failing="MCP selected integration tests"))
    steps = []
    precheck._mcp_checks({"mode": "affected", "base": "base"}, steps, tmp_path / "precheck.json")
    rows = {name: (cmd, env, cwd) for name, cmd, env, cwd in commands}
    assert "--base" in rows["MCP test selection"][0]
    cmd, env, cwd = rows["MCP selected package tests"]
    assert cmd[-1] == "tests" and "not xval and not slow" in cmd
    assert json.loads(env["RADIA_MCP_CI_SELECTION_JSON"]) == ["tests/test_ci_selection.py"]
    assert env["PYTEST_DISABLE_PLUGIN_AUTOLOAD"] == "1"
    assert cwd == ROOT / "packages/radia-mcp"
    assert "--confcutdir=tests/mcp_integration" in rows["MCP selected integration tests"][0]
    assert "MCP response tests" in rows and "MCP package test boundary" in rows
    assert any(step["returncode"] for step in steps)


def test_full_executes_prior_source_steps_and_full_mcp_plan(tmp_path, monkeypatch):
    commands = []
    monkeypatch.setattr(precheck, "_light_python", lambda *a: sys.executable)
    monkeypatch.setattr(precheck, "_run", fake_run(commands, selection=selection(package_tests=["tests"])))
    report = tmp_path / "report.json"
    assert precheck.main(["--full", "--report", str(report)]) == 0
    rows = {name: command for name, command, _, _ in commands}
    assert "--since" not in rows["fast contracts"] and "--full" in rows["MCP test selection"]
    assert "repository contracts" in rows and "CPU version" in rows
    assert json.loads(report.read_text())["qualified"]


def test_failed_selected_test_makes_main_fail(tmp_path, monkeypatch):
    plan = {"mode": "affected", "base": "base", "radia": False, "mcp": True}
    monkeypatch.setattr(precheck, "build_source_plan", lambda *a, **k: plan)
    monkeypatch.setattr(precheck, "_run", fake_run([], selection=selection(),
                                                 failing="MCP selected integration tests"))
    report = tmp_path / "report.json"
    assert precheck.main(["--report", str(report)]) == 1
    assert not json.loads(report.read_text())["qualified"]


@pytest.mark.parametrize("bad", ["../escape.py", "/tmp/escape.py", "tests/../../escape.py"])
def test_invalid_selected_paths_do_not_execute_tests(tmp_path, monkeypatch, bad):
    commands = []
    monkeypatch.setattr(precheck, "_run", fake_run(commands, selection=selection(package_tests=[bad])))
    with pytest.raises(ValueError):
        precheck._mcp_checks({"mode": "affected", "base": "base"}, [], tmp_path / "report.json")
    assert [row[0] for row in commands] == ["MCP test selection"]


def test_full_cannot_be_downgraded_to_a_partial_gate(tmp_path):
    with pytest.raises(SystemExit):
        precheck.main(["--full", "--skip-repo-contracts", "--report", str(tmp_path / "report.json")])


def test_caller_filters_cannot_silently_remove_required_tests(monkeypatch):
    for name in ("RADIA_MCP_CI_SELECTION_JSON", "RADIA_SOURCE_ONLY_CONTRACTS",
                 "RADIA_TESTS_ALLOW_PARTIAL", "PYTEST_ADDOPTS"):
        monkeypatch.setenv(name, "untrusted-caller-selection")
    env = precheck._clean_env()
    assert "RADIA_MCP_CI_SELECTION_JSON" not in env
    assert "RADIA_SOURCE_ONLY_CONTRACTS" not in env
    assert "PYTEST_ADDOPTS" not in env
    assert env["RADIA_TESTS_ALLOW_PARTIAL"] == "0"
