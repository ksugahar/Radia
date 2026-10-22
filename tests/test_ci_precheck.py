"""Precheck coverage must match CI and propagate real selected failures."""
import importlib.util
import json
import os
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


@pytest.mark.parametrize("parent_utf8", [None, "1"])
@pytest.mark.parametrize("exit_code", [0, 7])
def test_child_capture_is_independent_of_parent_encoding(
    tmp_path, parent_utf8, exit_code,
):
    child = tmp_path / "child.py"
    child.write_text(
        "import json, os, sys\n"
        "assert os.environ['PYTHONUTF8'] == '1'\n"
        "assert os.environ['PYTHONIOENCODING'] == 'utf-8'\n"
        "assert sys.stdout.encoding.lower() == 'utf-8'\n"
        "sys.stdout.buffer.write(bytes.fromhex('93fa967b8cea') + b'\\n')\n"
        "sys.stderr.buffer.write(bytes.fromhex('93fa967b8cea') + b'\\n')\n"
        "sys.stdout.flush(); sys.stderr.flush()\n"
        f"raise SystemExit({exit_code})\n",
        encoding="utf-8",
    )
    driver = tmp_path / "driver.py"
    driver.write_text(
        "import importlib.util, json, sys\n"
        f"spec=importlib.util.spec_from_file_location('precheck', {str(ROOT / 'tools/ci_precheck.py')!r})\n"
        "module=importlib.util.module_from_spec(spec); spec.loader.exec_module(module)\n"
        f"module.LOG_DIR=__import__('pathlib').Path({str(tmp_path / 'logs')!r})\n"
        "steps=[]\n"
        f"ok=module._run('legacy bytes', [sys.executable, {str(child)!r}], module._clean_env(), steps)\n"
        f"__import__('pathlib').Path({str(tmp_path / 'result.json')!r}).write_text(json.dumps(dict(ok=ok, steps=steps)), encoding='utf-8')\n",
        encoding="utf-8",
    )
    env = dict(os.environ)
    env.pop("PYTHONUTF8", None)
    env["PYTHONIOENCODING"] = "cp932:strict"
    if parent_utf8 is not None:
        env["PYTHONUTF8"] = parent_utf8
    result = subprocess.run([sys.executable, str(driver)], env=env,
                            capture_output=True, encoding="utf-8", errors="replace")
    assert result.returncode == 0, result.stdout + result.stderr
    recorded = json.loads((tmp_path / "result.json").read_text(encoding="utf-8"))
    assert recorded["ok"] is (exit_code == 0)
    assert recorded["steps"][0]["returncode"] == exit_code
    expected = bytes.fromhex("93fa967b8cea").decode("utf-8", errors="replace")
    log = Path(recorded["steps"][0]["log"]).read_text(encoding="utf-8")
    assert log == expected + "\n" + expected + "\n"


@pytest.mark.parametrize("exit_code", [0, 7])
def test_missing_capture_streams_keep_the_child_returncode(monkeypatch, exit_code):
    from types import SimpleNamespace

    def run(*args, **kwargs):
        assert kwargs["encoding"] == "utf-8"
        assert kwargs["errors"] == "replace"
        assert kwargs["env"]["PYTHONUTF8"] == "1"
        assert kwargs["env"]["PYTHONIOENCODING"] == "utf-8"
        return SimpleNamespace(returncode=exit_code, stdout=None, stderr=None)

    monkeypatch.setattr(precheck.subprocess, "run", run)
    steps = []
    assert precheck._run("empty streams", ["child"], {}, steps) is (exit_code == 0)
    assert steps[0]["returncode"] == exit_code
    assert steps[0]["tail"] == []


def git_diff(paths):
    def git(*args):
        if args[0] == "merge-base":
            return "a" * 40
        assert "--no-renames" in args and "-z" in args
        return "\0".join(paths) + ("\0" if paths else "")
    return git


@pytest.fixture(autouse=True)
def clean_candidate(monkeypatch):
    monkeypatch.setattr(precheck, "dirty_paths", lambda: [])


def workflow_commands(job):
    workflow = yaml.load((ROOT / ".github/workflows/policy-lint.yml").read_text(), Loader=yaml.BaseLoader)
    return {tuple(shlex.split(line.strip())[1:])
            for step in workflow["jobs"][job]["steps"]
            for line in step.get("run", "").splitlines()
            if line.strip().startswith("python ")
            and not line.strip().startswith("python -m pip install ")}


def test_policy_and_repository_commands_match_ci_in_both_directions():
    assert workflow_commands("policy-checks") == {args for _, args in precheck.POLICY_COMMANDS}
    assert workflow_commands("repository-contracts") == {args for _, args in precheck.REPOSITORY_COMMANDS}


def test_mcp_static_commands_match_workflow_in_both_directions():
    workflow = yaml.load((ROOT / ".github/workflows/radia-mcp-matrix.yml").read_text(), Loader=yaml.BaseLoader)
    names = ("Compile check (", "Publish-boundary lint", "Keep new gates out of")
    steps = [step for job in workflow["jobs"].values() for step in job["steps"]
             if step.get("name", "").startswith(names)]
    assert len(steps) == 3
    actual = set()
    for step in steps:
        for line in step["run"].splitlines():
            line = line.strip()
            if line.startswith(("python ", "python3 ")):
                parts = tuple(shlex.split(line)[1:])
                if "--base" in parts:
                    assert parts[-1] == "$CI_BASE_SHA"
                    actual.add((parts[0], "--base", "<base>"))
                else:
                    actual.add(parts)
    expected = {args for _, args in precheck.MCP_STATIC_COMMANDS}
    expected.add(("packages/radia-mcp/tools/check_ngsolve_architecture.py", "--base", "<base>"))
    assert actual == expected


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


def test_full_keeps_a_verified_base_and_empty_diff_keeps_source_baseline():
    assert precheck.build_source_plan("base", full=True, git=git_diff([]))["mcp"]
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
    assert [row[0] for row in commands] == [name for name, _ in (*precheck.POLICY_COMMANDS, *precheck.REPOSITORY_COMMANDS)]
    assert precheck.main(["--plan-only", "--report", str(report)]) == 0
    result = json.loads(report.read_text())
    assert result["status"] == "planned" and result["passed"] is None and not result["qualified"]


def test_invalid_base_records_failure_before_any_execution(tmp_path, monkeypatch, capsys):
    def unavailable(*args, **kwargs):
        raise subprocess.CalledProcessError(128, ["git", "merge-base"], stderr="fatal: missing comparison ref")
    monkeypatch.setattr(precheck, "build_source_plan", unavailable)
    monkeypatch.setattr(precheck, "_run", lambda *a, **k: pytest.fail("must not execute"))
    report = tmp_path / "report.json"
    assert precheck.main(["--report", str(report)]) == 1
    result = json.loads(report.read_text())
    assert not result["passed"]
    assert "missing comparison ref" in result["errors"][0]
    assert "missing comparison ref" in capsys.readouterr().out


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
    assert rows["MCP package test boundary"][0][1] == "packages/radia-mcp/tools/check_package_test_boundary.py"
    assert "--scratch-root" in rows["MCP package test boundary"][0]
    assert env["PYTHONPATH"] == str(ROOT / "packages/radia-mcp/src")
    assert "MCP collection dependencies" in rows
    assert any(step["returncode"] for step in steps)


def test_full_executes_prior_source_steps_and_full_mcp_plan(tmp_path, monkeypatch):
    # This is an orchestration test: the checked comparison is supplied by a
    # Git double, independent of checkout depth, remote refs and runner files.
    build_plan = precheck.build_source_plan
    monkeypatch.setattr(precheck, "build_source_plan",
                        lambda since, full=False: build_plan(since, full=full, git=git_diff([])))
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


@pytest.mark.parametrize("bad", ["../escape.py", "/tmp/escape.py", "tests/../../escape.py", r"tests\foo.py"])
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


def test_dirty_untracked_input_cannot_qualify(tmp_path, monkeypatch):
    monkeypatch.setattr(precheck, "dirty_paths", lambda: ["?? packages/radia-mcp/src/new.py"])
    monkeypatch.setattr(precheck, "_run", lambda *a, **k: pytest.fail("must not execute"))
    out = tmp_path / "report.json"
    assert precheck.main(["--report", str(out)]) == 1
    result = json.loads(out.read_text())
    assert result["dirty_paths"] and not result["qualified"]


def test_unexecuted_distribution_lanes_are_explicit():
    plan = precheck.build_source_plan("base", git=git_diff(["packages/cubit-mesh-export/src/cubit_mesh_export/__init__.py"]))
    assert "cubit-mesh-export" in plan["ci_only_lanes"]
    assert not plan["mcp"]


def test_missing_collection_dependencies_stop_before_pytest(tmp_path, monkeypatch):
    commands = []
    monkeypatch.setattr(precheck, "_run", fake_run(commands, selection=selection(),
                                                 failing="MCP collection dependencies"))
    precheck._mcp_checks({"mode": "affected", "base": "base"}, [], tmp_path / "report.json")
    assert "MCP collection dependencies" == commands[-1][0]
    assert not any("pytest" in cmd for _, cmd, _, _ in commands)


def test_partial_is_not_a_passing_qualification(tmp_path, monkeypatch):
    plan = precheck.build_source_plan("base", git=git_diff(["packages/cubit-mesh-export/README.md"]))
    monkeypatch.setattr(precheck, "build_source_plan", lambda *a, **k: plan)
    monkeypatch.setattr(precheck, "_run", fake_run([]))
    out = tmp_path / "report.json"
    assert precheck.main(["--skip-repo-contracts", "--report", str(out)]) == 0
    result = json.loads(out.read_text())
    assert result["passed"] is None and not result["qualified"] and result["status"] == "partial"


def test_actual_git_rename_and_delete_keep_old_and_new_sides(tmp_path):
    # Only the distribution diff promises both sides; MCP's selector is unchanged.
    def git(*args):
        return subprocess.check_output(["git", *args], cwd=tmp_path, text=True)
    git("init", "-q")
    git("config", "user.email", "test@example.invalid")
    git("config", "user.name", "Regression")
    (tmp_path / "src/radia").mkdir(parents=True)
    (tmp_path / "src/radia/old.py").write_text("x=1\n")
    (tmp_path / "src/radia/deleted.py").write_text("y=2\n")
    git("add", ".")
    git("commit", "-qm", "base")
    base = git("rev-parse", "HEAD").strip()
    (tmp_path / "docs").mkdir()
    git("mv", "src/radia/old.py", "docs/new.py")
    git("rm", "src/radia/deleted.py")
    git("commit", "-qm", "rename and delete")
    plan = precheck.build_source_plan(base, git=git)
    assert plan["changed_paths"] == ["docs/new.py", "src/radia/deleted.py", "src/radia/old.py"]


def test_boundary_helper_is_the_workflows_standalone_package_probe():
    workflow = (ROOT / ".github/workflows/radia-mcp-matrix.yml").read_text()
    assert 'sys.executable, "packages/radia-mcp/tools/check_package_test_boundary.py"' in workflow
    helper = (ROOT / "packages/radia-mcp/tools/check_package_test_boundary.py").read_text()
    assert 'env["RADIA_REPO_ROOT"] = str(root / "no-monorepo")' in helper
    assert 'TemporaryDirectory(' in helper


def test_default_report_sidecars_do_not_dirty_next_precheck(tmp_path, monkeypatch):
    # Use the real repository ignore rules and Git status, not a fake clean tree.
    subprocess.run(["git", "init", "-q", str(tmp_path)], check=True)
    (tmp_path / ".gitignore").write_bytes((ROOT / ".gitignore").read_bytes())
    subprocess.run(["git", "add", ".gitignore"], cwd=tmp_path, check=True)
    subprocess.run(["git", "-c", "user.email=test@example.invalid", "-c", "user.name=Test",
                    "commit", "-qm", "ignore policy"], cwd=tmp_path, check=True)
    monkeypatch.setattr(precheck, "ROOT", tmp_path)
    commands = []
    monkeypatch.setattr(precheck, "_run", fake_run(commands, selection=selection()))
    precheck._mcp_checks({"mode": "affected", "base": "base"}, [], tmp_path / "ci-precheck-report.json")
    artifacts = tmp_path / "ci-precheck-report-logs"
    # Collection diagnostics/helper scratch use the same ignored parent.
    (artifacts / "mcp-not-collected.json").write_text("{}")
    (artifacts / "boundary-scratch").mkdir()
    assert (artifacts / "mcp-selection.json").is_file()
    assert subprocess.check_output(["git", "status", "--porcelain",
                                    "--untracked-files=all"], cwd=tmp_path, text=True) == ""


def test_failed_partial_reports_false_not_unknown(tmp_path, monkeypatch):
    plan = precheck.build_source_plan("base", git=git_diff(["packages/cubit-mesh-export/README.md"]))
    monkeypatch.setattr(precheck, "build_source_plan", lambda *a, **k: plan)
    monkeypatch.setattr(precheck, "_run", fake_run([], failing="policy tests"))
    out = tmp_path / "report.json"
    assert precheck.main(["--skip-repo-contracts", "--report", str(out)]) == 1
    result = json.loads(out.read_text())
    assert result["passed"] is False and result["status"] == "failed"
