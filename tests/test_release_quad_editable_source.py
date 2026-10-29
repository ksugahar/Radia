from __future__ import annotations

import base64
import hashlib
import importlib.util
import os
import shutil
import subprocess
import sys
from pathlib import Path

import pytest


_TOOL = Path(__file__).resolve().parents[1] / "tools" / "release_quad.py"
_SPEC = importlib.util.spec_from_file_location("radia_release_quad_tool", _TOOL)
assert _SPEC is not None and _SPEC.loader is not None
release_quad = importlib.util.module_from_spec(_SPEC)
_SPEC.loader.exec_module(release_quad)
_GIT = shutil.which("git")
assert _GIT is not None, "Git is required by the release-quad contract tests"


@pytest.fixture(autouse=True)
def explicit_validation_runtime(monkeypatch):
    monkeypatch.setenv("RADIA_RELEASE_SSH_VALIDATION", "validation")
    monkeypatch.setenv("RADIA_RELEASE_PYTHON_VALIDATION", "C:/Runtime/python.exe")


def _git(repo, *args):
    return subprocess.run(
        [_GIT, "-C", str(repo), *args],
        check=True,
        capture_output=True,
        text=True,
    ).stdout.strip()


def _fake_wheel_site(root, *, shadow=False):
    """Create a tiny Radia wheel install with a valid hash-checked RECORD."""
    site = root / "site"
    package = site / "radia"
    dist_info = site / "radia-5.1.0.dist-info"
    package.mkdir(parents=True)
    dist_info.mkdir()
    files = {
        "radia/__init__.py": b"__version__ = '5.1.0'\n",
        "radia-5.1.0.dist-info/METADATA": (
            b"Metadata-Version: 2.1\nName: radia\nVersion: 5.1.0\n"
        ),
    }
    for relative, content in files.items():
        (site / relative).write_bytes(content)
    rows = []
    for relative, content in files.items():
        digest = base64.urlsafe_b64encode(hashlib.sha256(content).digest())
        rows.append(f"{relative},sha256={digest.rstrip(b'=').decode()},{len(content)}")
    rows.append("radia-5.1.0.dist-info/RECORD,,")
    (dist_info / "RECORD").write_text("\n".join(rows) + "\n", encoding="utf-8")
    paths = [site]
    if shadow:
        shadow_root = root / "shadow"
        shadow_package = shadow_root / "radia"
        shadow_package.mkdir(parents=True)
        (shadow_package / "__init__.py").write_bytes(files["radia/__init__.py"])
        paths.insert(0, shadow_root)
    env = os.environ.copy()
    env["PYTHONPATH"] = os.pathsep.join(str(path) for path in paths)
    return env


def test_release_quad_git_helper_trusts_only_its_active_worktree(monkeypatch):
    calls = []

    def fake_run(command, **kwargs):
        calls.append((command, kwargs))
        return subprocess.CompletedProcess(command, 0, stdout="", stderr="")

    monkeypatch.setattr(release_quad.subprocess, "run", fake_run)
    release_quad._git("status", "--short")

    command, kwargs = calls[0]
    assert command[:3] == [
        release_quad.GIT_EXE,
        "-c",
        f"safe.directory={release_quad.REPO.resolve().as_posix()}",
    ]
    assert command[3:] == [
        "-C", str(release_quad.REPO), "status", "--short"
    ]
    assert kwargs == {"capture_output": True, "text": True, "check": True}


def test_release_head_uses_version_tag_commit(monkeypatch):
    monkeypatch.setattr(release_quad, "_radia_version", lambda: "5.1.0")
    calls = []

    def fake_git(*args, **kwargs):
        calls.append((args, kwargs))
        return subprocess.CompletedProcess(
            args, 0, stdout="a" * 40 + "\n", stderr="")

    monkeypatch.setattr(release_quad, "_git", fake_git)

    assert release_quad._release_head() == "a" * 40
    assert calls == [(('rev-list', '-n', '1', 'v5.1.0'), {'check': False})]


@pytest.mark.parametrize("controller", ["build-machine", "validation-machine"])
@pytest.mark.parametrize("targets", ["validation,100", "all"])
def test_phase8_routes_explicit_hosts_independently_of_controller(monkeypatch, controller, targets):
    from argparse import Namespace

    calls = []
    monkeypatch.setattr(release_quad, "cmd_preflight", lambda _args: 0)
    monkeypatch.setattr(release_quad.platform, "node", lambda: controller)
    monkeypatch.setattr(release_quad, "_editable_repo_100", lambda: "release-source")
    monkeypatch.setattr(release_quad, "_deploy_pypi", lambda host, label, **kw: calls.append(host) or 0)
    monkeypatch.setattr(release_quad, "_deploy_editable_remote", lambda host, *args: calls.append((host, "development")) or 0)
    assert release_quad.cmd_phase8(Namespace(target=targets)) == 0
    expected = ["validation", release_quad.SSH_100, (release_quad.SSH_100, "development")]
    assert calls == expected


def test_default_editable_root_is_the_100_release_checkout(monkeypatch):
    monkeypatch.delenv(release_quad.EDITABLE_REPO_100_ENV, raising=False)
    monkeypatch.setattr(release_quad, "_radia_version", lambda: "5.1.0")
    monkeypatch.setattr(release_quad, "_release_commit", lambda: "c" * 40)
    assert release_quad._editable_repo_100() == (
        r"W:\00_CAE\Radia\release-quad\v5.1.0-ccccccccc")

def test_editable_runtime_gate_requires_exact_native_manifest():
    gate = release_quad.EDITABLE_RELEASE_VERIFY
    assert "release_native_payloads.json" in gate
    assert 'radia.release-native-payloads.v1' in gate
    assert '{"_radia_pybind.pyd", "axifem.pyd", "sparsesolv_ngsolve.pyd"}' in gate
    assert "native payload set differs" in gate
    # The retired native extension name is spelled in two pieces in the gate.
    assert 'name == "c" + "ln_core.pyd"' in gate
    assert '"locked-old" in lowered' in gate


def test_deployment_plan_uses_lab_wheel_and_100_development_editable(monkeypatch, capsys):
    monkeypatch.setattr(release_quad, "_release_commit", lambda: "c" * 40)
    monkeypatch.setattr(release_quad, "run", lambda *a, **k: pytest.fail("not a dry run"))
    monkeypatch.setattr(release_quad, "_pip_show", lambda *a: pytest.fail("runtime inference"))
    assert release_quad.cmd_deployment_plan(None) == 0
    plans = __import__("json").loads(capsys.readouterr().out)
    assert [plan["solver_runtime"] for plan in plans] == [
        "verified-wheel", "release-wheel+development-editable"]
    assert "development_source" not in plans[0]
    assert plans[1]["development_source"] == release_quad._editable_repo_100()
    assert all(plan["solver_action"] == "verify-and-install" for plan in plans)
    assert all(plan["independent_packages"] == "unchanged" for plan in plans)
    assert all(not plan["verified"] for plan in plans)
    assert all("mcp" not in key and "cubit" not in key
               for plan in plans for key in plan)

@pytest.mark.parametrize("rc", [0, 2, 4])
def test_lab_deploy_routes_wheel_to_lab_and_propagates_failure(monkeypatch, rc):
    calls = []
    monkeypatch.setattr(release_quad, "_deploy_pypi",
                        lambda host, label, **kw: calls.append((host, label, kw)) or rc)
    assert release_quad._deploy_validation() == rc
    assert calls == [("validation", "validation", {"python_path": "C:/Runtime/python.exe"})]

@pytest.mark.parametrize("drift", [0, 1])
def test_remote_deploy_changes_only_radia(monkeypatch, drift):
    monkeypatch.setattr(release_quad, "_release_head", lambda: "b" * 40)
    calls = []

    def verify(host, label, packages):
        assert host == "100" and packages == [("radia", "W:/release")]
        calls.append("verify")
        return drift

    monkeypatch.setattr(release_quad, "_verify_remote_editable", verify)
    def run(command, **_kwargs):
        calls.append(command)
        return subprocess.CompletedProcess(command, 0)
    monkeypatch.setattr(release_quad, "run", run)
    monkeypatch.setattr(
        release_quad, "_record_release_intent_remote",
        lambda host, label, repo: calls.append(("record", host, repo)) or 0)

    venv = release_quad.DEV_PYTHON_100
    assert release_quad._deploy_editable_remote(
        "100", "100", "W:/release", venv) == (4 if drift else 0)
    script = base64.b64decode(calls[0][-1]).decode("utf-16le")
    assert "pip install --no-deps" in script
    assert "pip uninstall" not in script
    assert "radia-mcp" not in script
    assert "cubit-mesh-export" not in script
    assert "Stop-Process" not in script
    assert "status --porcelain --untracked-files=no" in script
    assert "ngsolve.__version__" in script
    # Every interpreter call is the dedicated development venv, never PATH's python.
    assert f'& "{venv}" "-c"' in script
    assert f"& '{venv}' -m pip install --no-deps" in script
    assert script.index(f"Test-Path -LiteralPath '{venv}'") < script.index("pip install")
    assert '"python"' not in script and "\npython " not in script
    assert "netgen-mesher" in script
    assert script.index("rev-parse HEAD") < script.index("pip install")
    assert script.index("ngsolve.__version__") < script.index("pip install")
    if drift:
        assert calls[-1] == "verify"
        assert not any(isinstance(c, tuple) for c in calls)
    else:
        assert calls[-2] == "verify"
        assert calls[-1] == ("record", "100", "W:/release")


def test_done_checks_lab_wheel_and_100_editable_only(monkeypatch, tmp_path):
    from argparse import Namespace
    monkeypatch.setattr(release_quad, "_release_commit", lambda: "c" * 40)
    calls = []
    monkeypatch.setattr(release_quad, "cmd_preflight", lambda _args: 0)
    monkeypatch.setattr(release_quad, "cmd_temp_shadows", lambda _args: 0)
    monkeypatch.setattr(release_quad, "_release_head", lambda: "c" * 40)
    monkeypatch.setattr(release_quad, "_verify_local_release_source",
                        lambda root, sha: calls.append((root, sha)) or 0)
    monkeypatch.setattr(release_quad, "_verify_head_release_tag", lambda: calls.append("tag") or 0)
    monkeypatch.setattr(release_quad, "_verify_validation_wheel", lambda: calls.append("lab-wheel") or 0)
    monkeypatch.setattr(release_quad, "_verify_100_release_wheel",
                        lambda: calls.append("100-release") or 0)
    monkeypatch.setattr(release_quad, "_verify_100_editable",
                        lambda: calls.append("100-editable") or 0)
    monkeypatch.setattr(release_quad, "cmd_phase9", lambda _args: 4)
    assert release_quad.cmd_done(Namespace(simulink_package=None)) == 4
    assert calls == [(str(release_quad.REPO), "c" * 40), "tag", "lab-wheel",
                     "100-release", "100-editable"]

def test_cpu_release_does_not_require_matlab_candidate(monkeypatch, tmp_path):
    monkeypatch.setattr(release_quad, "_release_commit", lambda: "c" * 40)
    from argparse import Namespace

    monkeypatch.setattr(release_quad, "cmd_preflight", lambda _a: 0)
    monkeypatch.setattr(release_quad, "_verify_local_release_source", lambda *_a: 0)
    monkeypatch.setattr(release_quad, "cmd_temp_shadows", lambda _a: 0)
    monkeypatch.setattr(release_quad, "_verify_head_release_tag", lambda: 0)
    monkeypatch.setattr(release_quad, "_verify_validation_wheel", lambda: 0)
    monkeypatch.setattr(release_quad, "_verify_100_release_wheel", lambda: 0)
    monkeypatch.setattr(release_quad, "_verify_100_editable", lambda: 0)
    monkeypatch.setattr(release_quad, "cmd_phase9", lambda _a: 0)
    monkeypatch.setattr(release_quad, "_run_retired_standalone_pyside_guard", lambda: 0)
    monkeypatch.setattr(release_quad, "_check_main_synced", lambda **_k: 0)
    monkeypatch.setattr(release_quad, "_radia_version", lambda: "5.1.0")

    monkeypatch.setattr(release_quad, "_verify_final_pip_checks", lambda: 0)
    assert release_quad.cmd_done(Namespace(simulink_package=None)) == 0
    assert release_quad.cmd_done(Namespace(simulink_package="retired.zip")) == 4



@pytest.mark.parametrize("failure", [None, 0, 1, 2])
def test_done_pip_check_checks_all_declared_runtimes_over_ssh(monkeypatch, failure):
    calls = []

    def fake_run(command, **_kwargs):
        rc = int(len(calls) == failure)
        calls.append(command)
        return subprocess.CompletedProcess(command, rc)

    monkeypatch.setattr(release_quad, "run", fake_run)
    assert release_quad._verify_final_pip_checks() == (0 if failure is None else 4)
    assert [call[:2] for call in calls] == [
        ["ssh", "validation"], ["ssh", release_quad.SSH_100], ["ssh", release_quad.SSH_100]]
    remote = [base64.b64decode(call[-1]).decode("utf-16le") for call in calls]
    assert remote == [release_quad._runtime_guard_powershell("C:/Runtime/python.exe") + "& 'C:/Runtime/python.exe' -m pip check", f"& '{release_quad.DEV_PYTHON_100}' -m pip check", "python -m pip check"]


def test_100_editable_tooling_uses_the_dedicated_development_venv(monkeypatch):
    # 100号機's editable is a maintainer/student venv, separate from its release
    # runtime; the machine-wide python is never the editable target.
    assert release_quad.DEV_PYTHON_100 == (
        r"W:\00_CAE\Radia\environments\development\Scripts\python.exe")
    commands = []

    def run(command, source, timeout):
        commands.append(command)
        return subprocess.CompletedProcess(command, 0, "{}", "")

    monkeypatch.setattr(release_quad, "_run_script_with_file_output", run)
    release_quad._remote_editable_intent("100", ["--json", "verify"])
    release_quad._remote_editable_intent("mdx1", ["--json", "verify"])
    assert commands[0][:3] == ["ssh", "100", release_quad.DEV_PYTHON_100_PS]
    assert commands[1][:3] == ["ssh", "mdx1", "python"]


@pytest.mark.parametrize("wheel_rc", [0, 2, 3])
def test_100_deploy_updates_release_wheel_then_development_editable(monkeypatch, wheel_rc):
    # The machine-default Python is the student-facing release runtime (wheel);
    # the development venv holds the editable. A failed wheel stops both.
    called = []
    monkeypatch.setattr(release_quad, "_deploy_pypi",
                        lambda host, label, **kw: called.append(("wheel", host, label, kw)) or wheel_rc)
    monkeypatch.setattr(release_quad, "_deploy_editable_remote",
                        lambda *args: called.append(("editable", *args)) or 0)
    monkeypatch.setattr(release_quad, "_editable_repo_100", lambda: "W:/release")
    assert release_quad._deploy_100() == wheel_rc
    assert called[0] == ("wheel", "100", "build release runtime", {})
    if wheel_rc:
        assert len(called) == 1
    else:
        assert called[1] == ("editable", "100", "build", "W:/release",
                             release_quad.DEV_PYTHON_100)


@pytest.mark.skipif(shutil.which("pwsh") is None, reason="needs PowerShell 7 to execute the guard")
def test_missing_development_venv_exits_43_before_any_install(monkeypatch, tmp_path):
    # Execute the real PowerShell branch with a guaranteed-absent interpreter.
    absent = str(tmp_path / "absent venv" / "Scripts" / "python.exe")
    captured = {}
    monkeypatch.setattr(release_quad, "_release_head", lambda: "a" * 40)
    monkeypatch.setattr(release_quad, "run", lambda command, **_kw: captured.setdefault(
        "command", command) and subprocess.CompletedProcess(command, 1))
    assert release_quad._deploy_editable_remote("100", "100", str(tmp_path), absent) == 3
    result = subprocess.run(
        ["pwsh", "-NoProfile", "-NonInteractive", "-EncodedCommand", captured["command"][-1]],
        capture_output=True, text=True, timeout=60)
    assert result.returncode == release_quad.DEV_VENV_MISSING_EXIT == 43
    assert "Development venv interpreter is missing" in result.stderr
    assert "rev-parse" not in result.stdout + result.stderr



def test_restore_editable_is_a_tombstone_that_restores_nothing(monkeypatch, capsys):
    """The command stays only to refuse and name its replacement."""
    events = []
    monkeypatch.setattr(release_quad, "run",
                        lambda *a, **k: events.append(a) or None)
    assert release_quad.cmd_restore_editable(None) == 2
    out = capsys.readouterr().out
    assert "repoint" in out
    assert not events




def test_unc_normalization_covers_synthetic_canonical_and_release_worktrees(monkeypatch):
    monkeypatch.setattr(release_quad, "_PATH_ALIASES", (
        ("//192.0.2.10/share/repo/deep/", "x:/deep-root/"),
        ("//192.0.2.10/share/repo/", "x:/synthetic-root/"),
        ("//192.0.2.11/share/repo/", "x:/synthetic-root/"),
        ("w:/synthetic-repo/", "x:/synthetic-root/"),
    ))
    assert release_quad._norm_path(r"\\192.0.2.10\share\repo\deep\Candidate\\") == "x:/deep-root/candidate"
    assert release_quad._norm_path(r"\\192.0.2.10\share\repo\Release\src\radia") == "x:/synthetic-root/release/src/radia"
    assert release_quad._norm_path(r"\\192.0.2.11\SHARE\REPO\main\packages") == "x:/synthetic-root/main/packages"
    assert release_quad._norm_path(r"W:\synthetic-repo\main\\") == "x:/synthetic-root/main"
    assert release_quad._norm_path(r"\\unrelated-server\share\repo\main") == "//unrelated-server/share/repo/main"


def test_remote_default_aliases_match_local_normalization_contract():
    for alias, _ in release_quad._PATH_ALIASES[1:]:
        assert alias in release_quad.REMOTE_EDITABLE_VERIFY


def test_editable_probe_resolves_git_inside_the_target_process():
    probe = release_quad.CROSS_MACHINE_PROBE_100_EDITABLE

    assert 'git_exe = shutil.which("git")' in probe
    assert "[git_exe," in probe
    assert '"-c", "safe.directory=" + root' in probe
    assert "if result.returncode != 0:" in probe
    assert "git show failed for {relpath}" in probe
    assert "GIT_EXE" not in probe


def test_simulink_candidate_accepts_its_exact_tag_when_controller_is_newer(
        tmp_path, monkeypatch):
    repo = tmp_path / "release-controller"
    repo.mkdir()
    _git(repo, "init")
    _git(repo, "config", "user.name", "Radia Test")
    _git(repo, "config", "user.email", "radia-test@example.invalid")
    tracked = repo / "tracked.txt"
    tracked.write_text("candidate\n", encoding="ascii")
    _git(repo, "add", "tracked.txt")
    _git(repo, "commit", "-m", "candidate source")
    candidate = _git(repo, "rev-parse", "HEAD")
    _git(repo, "tag", "-a", "v4.95.75", "-m", "release", candidate)
    tracked.write_text("controller repair\n", encoding="ascii")
    _git(repo, "commit", "-am", "controller repair")
    head = _git(repo, "rev-parse", "HEAD")
    monkeypatch.setattr(release_quad, "REPO", repo)

    valid, message = release_quad._simulink_candidate_commit_is_release_anchored(
        {"commit": candidate, "radia_version": "4.95.75"}, head)

    assert valid
    assert "v4.95.75" in message

    valid, message = release_quad._simulink_candidate_commit_is_release_anchored(
        {"commit": candidate, "radia_version": "4.95.76"}, head)
    assert not valid
    assert "release tag" in message


def test_local_release_source_requires_exact_sha_and_tracked_clean(
        tmp_path, monkeypatch):
    repo = tmp_path / "release-source"
    repo.mkdir()
    _git(repo, "init")
    (repo / "tracked.txt").write_text("release\n", encoding="ascii")
    _git(repo, "add", "tracked.txt")
    _git(
        repo,
        "-c",
        "user.name=Radia Test",
        "-c",
        "user.email=radia-test@example.invalid",
        "commit",
        "-m",
        "release source",
    )
    head = _git(repo, "rev-parse", "HEAD")

    assert Path(release_quad.GIT_EXE).is_absolute()
    monkeypatch.setenv("PATH", "")
    assert release_quad._verify_local_release_source(str(repo), head) == 0
    assert release_quad._verify_local_release_source(str(repo), "0" * 40) == 4

    (repo / "tracked.txt").write_text("parallel WIP\n", encoding="ascii")
    assert release_quad._verify_local_release_source(str(repo), head) == 4


def _build_release_repo(root):
    """A tagged v5.2.1 release commit followed by a tooling repair commit."""
    repo = root / "controller"
    repo.mkdir()
    _git(repo, "init")
    _git(repo, "config", "user.name", "Radia Test")
    _git(repo, "config", "user.email", "radia-test@example.invalid")
    (repo / "pyproject.toml").write_text('# Release metadata — 日本語\n[project]\nname = "radia"\nversion = "5.2.1"\n',
                                         encoding="utf-8")
    (repo / "tool.py").write_text("release\n", encoding="ascii")
    _git(repo, "add", "pyproject.toml", "tool.py")
    _git(repo, "commit", "-m", "release")
    release = _git(repo, "rev-parse", "HEAD")
    _git(repo, "tag", "-a", "v5.2.1", "-m", "release", release)
    (repo / "tool.py").write_text("repaired release tooling\n", encoding="ascii")
    _git(repo, "commit", "-am", "tooling repair")
    return repo, release


@pytest.fixture(scope="module")
def release_template(tmp_path_factory):
    """Build the release controller and a clone of it once per module.

    Every test works on its own copy (a directory copy keeps the commits,
    tag and SHAs), so the per-test Git process cost is only what the test
    itself exercises.
    """
    root = tmp_path_factory.mktemp("release-template")
    repo, release = _build_release_repo(root)
    source = root / "release-source"
    _git(root, "clone", "-q", str(repo), str(source))
    return repo, source, release


def _release_repo(tmp_path, template):
    repo, _source, release = template
    copy = tmp_path / "controller"
    shutil.copytree(repo, copy)
    return copy, release


def _release_source(tmp_path, template):
    _repo, source, _release = template
    copy = tmp_path / "release-source"
    shutil.copytree(source, copy)
    return copy


def test_later_tooling_controller_must_descend_cleanly_from_the_release(
        tmp_path, release_template, monkeypatch):
    repo, release = _release_repo(tmp_path, release_template)
    monkeypatch.setattr(subprocess, "_text_encoding", lambda: "cp932")
    assert release_quad._verify_release_controller(repo, release, "5.2.1") == 0
    # Same version declared at the release commit is part of the identity.
    assert release_quad._verify_release_controller(repo, release, "5.2.2") == 4
    # Tracked changes in the controller are refused.
    (repo / "tool.py").write_text("unreviewed\n", encoding="ascii")
    assert release_quad._verify_release_controller(repo, release, "5.2.1") == 4
    _git(repo, "checkout", "--", "tool.py")
    # A controller that does not contain the release commit is refused.
    _git(repo, "checkout", "--orphan", "unrelated")
    _git(repo, "commit", "-m", "unrelated history")
    assert release_quad._verify_release_controller(repo, release, "5.2.1") == 4


def _done_until_source_gates(monkeypatch, controller, release):
    monkeypatch.setattr(release_quad, "REPO", controller)
    monkeypatch.setattr(release_quad, "cmd_preflight", lambda _args: 0)
    monkeypatch.setattr(release_quad, "_release_head", lambda: release)
    monkeypatch.setattr(release_quad, "_radia_version", lambda: "5.2.1")
    passed = AssertionError("source gates passed")

    def next_gate():
        raise passed

    monkeypatch.setattr(release_quad, "_verify_head_release_tag", next_gate)
    return passed


def test_done_accepts_a_later_controller_with_a_separate_exact_tag_source(
        tmp_path, monkeypatch, release_template):
    from argparse import Namespace
    controller, release = _release_repo(tmp_path, release_template)
    source = _release_source(tmp_path, release_template)
    _git(source, "checkout", "-q", "--detach", "v5.2.1")
    passed = _done_until_source_gates(monkeypatch, controller, release)
    with pytest.raises(AssertionError) as raised:
        release_quad.cmd_done(Namespace(simulink_package=None, release_source=str(source)))
    assert raised.value is passed
    # Without a separate source the later controller is itself the source: refused.
    assert release_quad.cmd_done(Namespace(simulink_package=None, release_source=None)) == 4


@pytest.mark.parametrize("case", ["stale-source", "dirty-source", "dirty-controller",
                                  "unrelated-controller"])
def test_done_rejects_a_bad_source_or_controller(tmp_path, monkeypatch, case, release_template):
    from argparse import Namespace
    controller, release = _release_repo(tmp_path, release_template)
    source = _release_source(tmp_path, release_template)
    if case != "stale-source":
        _git(source, "checkout", "-q", "--detach", "v5.2.1")
    if case == "dirty-source":
        (source / "tool.py").write_text("parallel WIP\n", encoding="ascii")
    elif case == "dirty-controller":
        (controller / "tool.py").write_text("unreviewed\n", encoding="ascii")
    elif case == "unrelated-controller":
        _git(controller, "checkout", "-q", "--orphan", "unrelated")
        _git(controller, "commit", "-q", "-m", "unrelated history")
    _done_until_source_gates(monkeypatch, controller, release)
    assert release_quad.cmd_done(
        Namespace(simulink_package=None, release_source=str(source))) == 4


@pytest.mark.skipif(shutil.which("pwsh") is None, reason="needs PowerShell 7 to execute the guard")
@pytest.mark.parametrize("case,code", [("wrong-sha", 41), ("dirty", 42)])
def test_editable_deploy_source_refusals_exit_with_their_codes(
        tmp_path, monkeypatch, case, code, release_template):
    import sys
    controller, release = _release_repo(tmp_path, release_template)
    expected = release if case == "dirty" else "0" * 40
    _git(controller, "checkout", "-q", "--detach", release)
    if case == "dirty":
        (controller / "tool.py").write_text("parallel WIP\n", encoding="ascii")
    captured = {}
    monkeypatch.setattr(release_quad, "_release_head", lambda: expected)
    monkeypatch.setattr(release_quad, "run", lambda command, **_kw: captured.setdefault(
        "command", command) and subprocess.CompletedProcess(command, 1))
    assert release_quad._deploy_editable_remote(
        "100", "100", str(controller), sys.executable) == 3
    result = subprocess.run(
        ["pwsh", "-NoProfile", "-NonInteractive", "-EncodedCommand", captured["command"][-1]],
        capture_output=True, text=True, timeout=60)
    assert result.returncode == code
    assert ("SHA mismatch" if code == 41 else "tracked changes") in result.stderr
    assert "pip install" not in result.stdout


def test_release_tag_gate_requires_declared_version_at_exact_head(monkeypatch):
    monkeypatch.setattr(
        release_quad, "_read_repo_versions", lambda: {"radia": "4.95.80"}
    )
    monkeypatch.setattr(release_quad, "_release_head", lambda: "a" * 40)
    monkeypatch.setattr(release_quad, "_release_tag_commit", lambda _version: "a" * 40)
    assert release_quad._verify_head_release_tag() == 0

    monkeypatch.setattr(release_quad, "_release_tag_commit", lambda _version: None)
    assert release_quad._verify_head_release_tag() == 4

    monkeypatch.setattr(release_quad, "_release_tag_commit", lambda _version: "b" * 40)
    assert release_quad._verify_head_release_tag() == 4


def test_lab_wheel_deploy_stops_before_ssh_when_unpublished(monkeypatch):
    monkeypatch.setattr(release_quad, "_read_repo_versions", lambda *_: {"radia": "5.2.0"})
    monkeypatch.setattr(release_quad, "_check_pypi_propagation", lambda *_: 2)
    monkeypatch.setattr(release_quad, "run", lambda *_a, **_k: pytest.fail("No SSH before wheel availability"))
    assert release_quad._deploy_validation() == 2


def test_lab_wheel_deploy_has_guard_no_editable_and_pip_check(monkeypatch):
    calls = []
    monkeypatch.setattr(release_quad, "_read_repo_versions", lambda *_: {"radia": "5.1.0"})
    monkeypatch.setattr(release_quad, "_check_pypi_propagation", lambda *_: 0)
    monkeypatch.setattr(release_quad, "run", lambda command, **kw: calls.append(command) or subprocess.CompletedProcess(command, 0))
    assert release_quad._deploy_validation() == 0
    assert calls[0][:2] == ["ssh", "validation"]
    script = base64.b64decode(calls[0][-1]).decode("utf-16le")
    assert script.index("base64") < script.index("pip install")
    assert "--only-binary=:all:" in script
    assert "--no-deps" in script
    assert " -e " not in script
    assert "direct_url.json" in script and "editable" in script
    assert "pip check" in script
    assert "radia-mcp" not in script and "cubit-mesh-export" not in script


def test_remote_deploy_checks_exact_source_before_install(monkeypatch):
    expected_sha = "a" * 40
    captured = {}
    monkeypatch.setattr(release_quad, "_release_head", lambda: expected_sha)

    def capture_run(command, **_kwargs):
        captured["command"] = command
        return subprocess.CompletedProcess(command, 0)

    monkeypatch.setattr(release_quad, "run", capture_run)
    monkeypatch.setattr(release_quad, "_verify_remote_editable", lambda *_a: 0)
    monkeypatch.setattr(release_quad, "_record_release_intent_remote", lambda *_a: 0)

    assert release_quad._deploy_editable_remote(
        "release-host", "release host", r"W:\Radia\release-source",
        release_quad.DEV_PYTHON_100
    ) == 0
    script = base64.b64decode(captured["command"][-1]).decode("utf-16le")
    assert expected_sha in script
    assert 'safe.directory=W:/Radia/release-source' in script
    assert "status --porcelain --untracked-files=no" in script
    assert "ngsolve.__version__" in script
    assert "netgen-mesher" in script
    assert script.index("rev-parse HEAD") < script.index("pip install --no-deps")
    assert script.index("ngsolve.__version__") < script.index("pip install --no-deps")
    assert "pip uninstall" not in script
    assert "radia-mcp" not in script
    assert "cubit-mesh-export" not in script


def test_done_keeps_lab_wheel_and_100_editable_checks_after_all_gates(monkeypatch):
    calls = []
    monkeypatch.setattr(release_quad, "cmd_temp_shadows", lambda _args: calls.append("shadows") or 0)
    monkeypatch.setattr(release_quad, "cmd_preflight", lambda _args: calls.append("preflight") or 0)
    monkeypatch.setattr(release_quad, "_release_head", lambda: "a" * 40)
    monkeypatch.setattr(release_quad, "_verify_local_release_source",
                        lambda _repo, _sha: calls.append("source") or 0)
    monkeypatch.setattr(release_quad, "_verify_head_release_tag", lambda: calls.append("tag") or 0)
    monkeypatch.setattr(release_quad, "_verify_validation_wheel", lambda *_args: calls.append("validation") or 0)
    monkeypatch.setattr(release_quad, "_verify_100_release_wheel",
                        lambda *_args: calls.append("100-release") or 0)
    monkeypatch.setattr(release_quad, "_verify_100_editable", lambda *_args: calls.append("100") or 0)
    monkeypatch.setattr(release_quad, "cmd_phase9", lambda _args: calls.append("phase9") or 0)
    monkeypatch.setattr(release_quad, "_run_retired_standalone_pyside_guard",
                        lambda: calls.append("guard") or 0)
    monkeypatch.setattr(release_quad, "_check_main_synced", lambda **_kwargs: calls.append("main") or 0)
    monkeypatch.setattr(release_quad, "_verify_simulink_candidate_state",
                        lambda package: calls.append(("simulink", package)) or 0)
    monkeypatch.setattr(release_quad, "_verify_final_pip_checks", lambda: calls.append("pip-check") or 0)
    args = type("Args", (), {"simulink_package": None})()
    assert release_quad.cmd_done(args) == 0
    assert calls == ["preflight", "source", "tag", "validation", "100-release", "100",
                     "phase9", "guard", "main", "pip-check"]

def test_done_stops_before_machine_checks_when_active_source_is_stale(monkeypatch):
    calls = []
    monkeypatch.setattr(release_quad, "cmd_preflight", lambda _args: calls.append("preflight") or 0)
    monkeypatch.setattr(release_quad, "_release_head", lambda: "a" * 40)
    monkeypatch.setattr(release_quad, "_verify_local_release_source",
                        lambda _root, _sha: calls.append("source") or 4)
    def unexpected_machine_check(*_args, **_kwargs):
        raise AssertionError("machine checks must not run from a stale source")
    monkeypatch.setattr(release_quad, "_verify_validation_wheel", unexpected_machine_check)
    monkeypatch.setattr(release_quad, "_verify_100_release_wheel", unexpected_machine_check)
    monkeypatch.setattr(release_quad, "_verify_100_editable", unexpected_machine_check)
    args = type("Args", (), {"simulink_package": None})()
    assert release_quad.cmd_done(args) == 4
    assert calls == ["preflight", "source"]

def test_ci_check_runs_use_latest_attempt_for_each_name(monkeypatch):
    runs = {
        "check_runs": [
            {
                "id": 1,
                "name": "fast-contracts",
                "status": "completed",
                "conclusion": "cancelled",
            },
            {
                "id": 2,
                "name": "fast-contracts",
                "status": "completed",
                "conclusion": "success",
            },
            {
                "id": 3,
                "name": "policy-checks",
                "status": "completed",
                "conclusion": "success",
            },
        ]
    }
    monkeypatch.setattr(release_quad, "_git_repo_owner_name", lambda: "owner/repo")
    monkeypatch.setattr(release_quad, "gh_get", lambda _path: (runs, {}))

    ok, message = release_quad._check_github_hosted_workflows(
        "a" * 40, required_names=None, timeout_sec=1, poll_sec=0
    )

    assert ok is True
    assert "2 latest SHA-bound check-runs GREEN" in message


def _green(name, run_id):
    return {"id": run_id, "name": name,
            "status": "completed", "conclusion": "success"}


def test_absent_release_check_is_not_passing_evidence(monkeypatch):
    """A SHA the release workflow never ran on must not verify.

    Radia Native Release triggers only on a v* tag or a manual dispatch, so a
    commit can be green on every check that did run while the release build
    has never seen it. Treating that silence as success is what spent tags
    4.95.86..4.95.89.
    """
    runs = {"check_runs": [_green("policy-checks", 1), _green("mcp-matrix", 2)]}
    monkeypatch.setattr(release_quad, "_git_repo_owner_name", lambda: "owner/repo")
    monkeypatch.setattr(release_quad, "gh_get", lambda _path: (runs, {}))

    ok, message = release_quad._check_github_hosted_workflows(
        "a" * 40, required_names=None,
        require_present=[release_quad.RELEASE_CHECK_RUN],
        timeout_sec=1, poll_sec=0, registration_grace_sec=0,
    )

    assert ok is False
    assert release_quad.RELEASE_CHECK_RUN in message


def test_present_release_check_verifies_with_the_rest(monkeypatch):
    """require_present adds a requirement; it does not narrow the check set."""
    runs = {"check_runs": [
        _green("policy-checks", 1),
        _green(release_quad.RELEASE_CHECK_RUN, 2),
    ]}
    monkeypatch.setattr(release_quad, "_git_repo_owner_name", lambda: "owner/repo")
    monkeypatch.setattr(release_quad, "gh_get", lambda _path: (runs, {}))

    ok, message = release_quad._check_github_hosted_workflows(
        "a" * 40, required_names=None,
        require_present=[release_quad.RELEASE_CHECK_RUN],
        timeout_sec=1, poll_sec=0,
    )

    assert ok is True
    assert "2 latest SHA-bound check-runs GREEN" in message


def test_require_present_still_fails_on_an_unrelated_red(monkeypatch):
    """Requiring the release build by name must not excuse any other check."""
    runs = {"check_runs": [
        _green(release_quad.RELEASE_CHECK_RUN, 1),
        {"id": 2, "name": "policy-checks",
         "status": "completed", "conclusion": "failure"},
    ]}
    monkeypatch.setattr(release_quad, "_git_repo_owner_name", lambda: "owner/repo")
    monkeypatch.setattr(release_quad, "gh_get", lambda _path: (runs, {}))

    ok, _message = release_quad._check_github_hosted_workflows(
        "a" * 40, required_names=None,
        require_present=[release_quad.RELEASE_CHECK_RUN],
        timeout_sec=1, poll_sec=0,
    )

    assert ok is False


def test_release_check_run_name_matches_the_workflow_job():
    """RELEASE_CHECK_RUN must name a job that build-test.yml actually defines."""
    workflow = (Path(__file__).resolve().parents[1]
                / ".github" / "workflows" / "build-test.yml").read_text(encoding="utf-8")
    assert f"\n    name: {release_quad.RELEASE_CHECK_RUN}\n" in workflow
    assert release_quad.RELEASE_CHECK_RUN != "build-test"


def test_other_distribution_build_does_not_certify_native_release(monkeypatch):
    runs = {"check_runs": [_green("build-test", 1), _green("policy-checks", 2)]}
    monkeypatch.setattr(release_quad, "_git_repo_owner_name", lambda: "owner/repo")
    monkeypatch.setattr(release_quad, "gh_get", lambda _path: (runs, {}))
    ok, message = release_quad._check_github_hosted_workflows(
        "a" * 40, require_present=[release_quad.RELEASE_CHECK_RUN],
        timeout_sec=1, poll_sec=0, registration_grace_sec=0,
    )
    assert not ok
    assert release_quad.RELEASE_CHECK_RUN in message


@pytest.mark.parametrize("conclusion", ["skipped", "neutral", "failure", "cancelled"])
@pytest.mark.parametrize("required_names", [None, ["policy-checks"]])
def test_native_evidence_must_succeed_even_when_other_checks_are_selected(
    monkeypatch, conclusion, required_names
):
    native = _green(release_quad.RELEASE_CHECK_RUN, 2)
    native["conclusion"] = conclusion
    runs = {"check_runs": [_green("policy-checks", 1), native]}
    monkeypatch.setattr(release_quad, "_git_repo_owner_name", lambda: "owner/repo")
    monkeypatch.setattr(release_quad, "gh_get", lambda _path: (runs, {}))
    ok, message = release_quad._check_github_hosted_workflows(
        "a" * 40, required_names=required_names,
        require_present=[release_quad.RELEASE_CHECK_RUN],
        timeout_sec=1, poll_sec=0, registration_grace_sec=0,
    )
    assert not ok
    assert f"{release_quad.RELEASE_CHECK_RUN}: {conclusion}" in message


@pytest.mark.parametrize("workflow_path", sorted(release_quad.INDEPENDENT_DISTRIBUTION_WORKFLOWS))
@pytest.mark.parametrize("status,conclusion", [("queued", None), ("completed", "failure")])
def test_cpu_release_does_not_wait_for_independent_distribution_ci(
    monkeypatch, workflow_path, status, conclusion
):
    sha = "a" * 40
    native = _green(release_quad.RELEASE_CHECK_RUN, 1)
    cpu = _green("fast-contracts", 2)
    independent = dict(id=3, name="fast-contracts", status=status, conclusion=conclusion,
        details_url="https://github.com/owner/repo/actions/runs/50/job/60")
    def github(path):
        if "/check-runs?" in path:
            return {"check_runs": [native, cpu, independent]}, {}
        assert path == "repos/owner/repo/actions/runs/50"
        return dict(head_sha=sha, path=workflow_path), {}
    monkeypatch.setattr(release_quad, "_git_repo_owner_name", lambda: "owner/repo")
    monkeypatch.setattr(release_quad, "gh_get", github)
    monkeypatch.setattr(release_quad, "_release_head", lambda: sha)
    assert release_quad.cmd_ci_verify(None) == 0


@pytest.mark.parametrize("workflow_path", [".github/workflows/radia.yml", "unknown-shared-check.yml"])
def test_cpu_release_still_rejects_cpu_and_unclassified_failures(monkeypatch, workflow_path):
    sha = "a" * 40
    failed = dict(id=2, name="shared-contracts", status="completed", conclusion="failure",
        details_url="https://github.com/owner/repo/actions/runs/50/job/60")
    def github(path):
        if "/check-runs?" in path:
            return {"check_runs": [_green(release_quad.RELEASE_CHECK_RUN, 1), failed]}, {}
        return dict(head_sha=sha, path=workflow_path), {}
    monkeypatch.setattr(release_quad, "_git_repo_owner_name", lambda: "owner/repo")
    monkeypatch.setattr(release_quad, "gh_get", github)
    ok, message = release_quad._check_github_hosted_workflows(sha,
        require_present=[release_quad.RELEASE_CHECK_RUN],
        excluded_workflows=release_quad.INDEPENDENT_DISTRIBUTION_WORKFLOWS,
        timeout_sec=0, poll_sec=0)
    assert not ok and "shared-contracts: failure" in message


@pytest.mark.parametrize("ownership", ["lookup-error", "wrong-sha", "missing-path"])
def test_cpu_release_cannot_ignore_unverified_workflow_ownership(monkeypatch, ownership):
    sha = "a" * 40
    gpu = dict(id=2, name="cuda-mex", status="completed", conclusion="failure",
        details_url="https://github.com/owner/repo/actions/runs/50/job/60")
    def github(path):
        if "/check-runs?" in path:
            return {"check_runs": [_green(release_quad.RELEASE_CHECK_RUN, 1), gpu]}, {}
        if ownership == "lookup-error":
            raise RuntimeError("ownership API unavailable")
        return dict(head_sha="b" * 40 if ownership == "wrong-sha" else sha,
            path=None if ownership == "missing-path" else ".github/workflows/radia-simulink.yml"), {}
    monkeypatch.setattr(release_quad, "_git_repo_owner_name", lambda: "owner/repo")
    monkeypatch.setattr(release_quad, "gh_get", github)
    ok, _ = release_quad._check_github_hosted_workflows(sha,
        require_present=[release_quad.RELEASE_CHECK_RUN],
        excluded_workflows=release_quad.INDEPENDENT_DISTRIBUTION_WORKFLOWS,
        timeout_sec=0, poll_sec=0)
    assert not ok


def test_gpu_named_native_check_cannot_replace_cpu_native_acceptance(monkeypatch):
    sha = "a" * 40
    gpu = _green(release_quad.RELEASE_CHECK_RUN, 1)
    gpu["details_url"] = "https://github.com/owner/repo/actions/runs/50/job/60"
    def github(path):
        if "/check-runs?" in path:
            return {"check_runs": [gpu]}, {}
        return dict(head_sha=sha, path=".github/workflows/radia-simulink.yml"), {}
    monkeypatch.setattr(release_quad, "_git_repo_owner_name", lambda: "owner/repo")
    monkeypatch.setattr(release_quad, "gh_get", github)
    ok, message = release_quad._check_github_hosted_workflows(sha,
        require_present=[release_quad.RELEASE_CHECK_RUN],
        excluded_workflows=release_quad.INDEPENDENT_DISTRIBUTION_WORKFLOWS,
        timeout_sec=0, poll_sec=0, registration_grace_sec=0)
    assert not ok and release_quad.RELEASE_CHECK_RUN in message


@pytest.mark.parametrize("host", ["100"])
def test_explicit_editable_root_does_not_resolve_controller_version_tag(monkeypatch, host):
    env = getattr(release_quad, "EDITABLE_REPO_" + host.upper() + "_ENV")
    monkeypatch.setenv(env, "C:/release-quad/v5.1.0-ccccccccc/")
    monkeypatch.setattr(release_quad, "_release_commit",
                        lambda: pytest.fail("explicit root must not resolve a controller tag"))
    assert getattr(release_quad, "_editable_repo_" + host)() == "C:/release-quad/v5.1.0-ccccccccc"


@pytest.mark.parametrize("verifier,host,label", [
    ("_verify_validation_wheel", "validation", "validation"),
    ("_verify_100_release_wheel", "100", "build release runtime"),
])
@pytest.mark.parametrize("rc,expected", [(0, 0), (1, 4)])
def test_wheel_verifiers_use_remote_record_gate(monkeypatch, verifier, host, label, rc, expected):
    calls = []
    monkeypatch.setattr(release_quad, "_radia_version", lambda: "5.1.0")
    monkeypatch.setattr(release_quad, "run", lambda cmd, **kw: calls.append(cmd) or subprocess.CompletedProcess(cmd, rc))
    assert getattr(release_quad, verifier)() == expected
    assert calls[0][:2] == ["ssh", host]
    script = base64.b64decode(calls[0][-1]).decode("utf-16le")
    # The release runtime is the machine-default interpreter, not a venv.
    assert ("& 'C:/Runtime/python.exe' -c" in script) if verifier == '_verify_validation_wheel' else script.startswith('python -c')
    payload = script.rsplit("b64decode('", 1)[1].split("'", 1)[0]
    code = base64.b64decode(payload).decode()
    assert "d.files" in code and "actual == f.hash.value" in code
    assert "Shadowed Radia import" in code and f"{label} must use a wheel" in code
    assert "'5.1.0'" in code


@pytest.mark.parametrize(
    "expected_version,shadow,wrong_root,expected_rc,expected_message",
    [
        ("5.1.0", False, False, 0, "validation wheel verified"),
        ("5.1.1", False, False, 1, "Radia version mismatch"),
        ("5.1.0", True, False, 1, "Shadowed Radia import"),
        ("5.1.0", False, True, 1, "Radia wheel must be under site-packages"),
    ],
)
def test_wheel_runtime_probe_checks_installed_version_and_import_origin(
    tmp_path, expected_version, shadow, wrong_root, expected_rc, expected_message
):
    env = _fake_wheel_site(tmp_path, shadow=shadow)
    code = release_quad._wheel_runtime_probe_code(expected_version, "validation")
    # Model the interpreter's install scheme without installing a test package
    # into the real environment; the production probe itself runs unchanged.
    install_root = str(tmp_path / ("other-site" if wrong_root else "site"))
    setup = (
        "import sysconfig\n"
        f"sysconfig.get_paths = lambda: {{'purelib': {install_root!r}, 'platlib': {install_root!r}}}\n"
    )

    result = subprocess.run(
        [sys.executable, "-c", setup + code],
        capture_output=True,
        text=True,
        env=env,
        timeout=30,
    )

    assert result.returncode == expected_rc
    assert expected_message in result.stdout + result.stderr


@pytest.mark.parametrize("status,conclusion", [("queued", None), ("completed", "failure")])
def test_cpu_release_excludes_a100_nightly(monkeypatch, status, conclusion):
    sha = "a" * 40
    gpu = dict(id=2, name="validate", status=status, conclusion=conclusion,
        details_url="https://github.com/owner/repo/actions/runs/50/job/60")
    def github(path):
        if "/check-runs?" in path:
            return {"check_runs": [_green(release_quad.RELEASE_CHECK_RUN, 1), gpu]}, {}
        return dict(head_sha=sha, path=".github/workflows/gpu-validation.yml"), {}
    monkeypatch.setattr(release_quad, "_git_repo_owner_name", lambda: "owner/repo")
    monkeypatch.setattr(release_quad, "gh_get", github)
    ok, message = release_quad._check_github_hosted_workflows(sha,
        require_present=[release_quad.RELEASE_CHECK_RUN],
        excluded_workflows=release_quad.INDEPENDENT_DISTRIBUTION_WORKFLOWS,
        timeout_sec=0, poll_sec=0)
    assert ok, message


@pytest.mark.parametrize("missing", ["RADIA_RELEASE_SSH_VALIDATION", "RADIA_RELEASE_PYTHON_VALIDATION"])
def test_validation_config_fails_before_any_phase8_action(monkeypatch, missing):
    """Guards: no fallback or partial deployment. Assumes: remote actions are observed mocks."""
    from argparse import Namespace
    monkeypatch.delenv(missing)
    monkeypatch.setattr(release_quad, "run", lambda *a, **kw: pytest.fail("No remote action"))
    monkeypatch.setattr(release_quad, "cmd_preflight", lambda *a: pytest.fail("No preflight before config"))
    with pytest.raises(ValueError, match=missing):
        release_quad.cmd_phase8(Namespace(target="100,validation"))


@pytest.mark.parametrize("name,value", [
    ("RADIA_RELEASE_SSH_VALIDATION", "-oProxyCommand=bad"),
    ("RADIA_RELEASE_SSH_VALIDATION", "target extra"),
    ("RADIA_RELEASE_PYTHON_VALIDATION", "python"),
    ("RADIA_RELEASE_PYTHON_VALIDATION", "C:relative/python.exe"),
    ("RADIA_RELEASE_PYTHON_VALIDATION", "C:/WindowsApps/python.exe"),
    ("RADIA_RELEASE_PYTHON_VALIDATION", 'C:/x/python.exe";bad'),
])
def test_validation_rejects_ambiguous_runtime_configuration(monkeypatch,name,value):
    """Guards: explicit executable/alias only. Assumes: Windows target semantics."""
    monkeypatch.setenv(name,value)
    with pytest.raises(ValueError,match=name):release_quad._solver_target("validation")


def test_all_validation_consumers_resolve_same_explicit_runtime(monkeypatch):
    """Guards: deploy/check/pip/phase9 share configuration. Assumes: no SSH executed."""
    import ntpath
    host="consumer-test";python="C:/Program Files/Owner's Python/python.exe"
    monkeypatch.setenv("RADIA_RELEASE_SSH_VALIDATION",host)
    monkeypatch.setenv("RADIA_RELEASE_PYTHON_VALIDATION",python)
    assert release_quad._solver_target("validation")==("validation",host,python)
    calls=[]
    monkeypatch.setattr(release_quad,"_deploy_pypi",lambda h,l,**kw:calls.append((h,kw)) or 0)
    assert release_quad._deploy_validation()==0
    assert calls==[(host,{"python_path":python})]
    commands=[]
    monkeypatch.setattr(release_quad,"_radia_version",lambda:"5.4.0")
    monkeypatch.setattr(release_quad,"run",lambda cmd,**kw:commands.append(cmd) or subprocess.CompletedProcess(cmd,0))
    assert release_quad._verify_validation_wheel()==0
    assert release_quad._verify_final_pip_checks()==0
    assert commands[0][1]==host and commands[1][1]==host
    literal="'C:/Program Files/Owner''s Python/python.exe'"
    for cmd in commands[:2]:
        script=base64.b64decode(cmd[-1]).decode("utf-16le")
        assert "Test-Path -LiteralPath "+literal in script,script
        assert "& "+literal+" -c" in script,script
        assert "Selected Python is missing" in script,script
        token=script.split("b64decode('",1)[1].split("'",1)[0]
        code=base64.b64decode(token).decode()
        assert "Selected Python identity mismatch" in code
        namespace={"__name__":"runtime_identity_probe"}
        # Independent oracle: production probe rejects an executable different from selection.
        with pytest.raises(AssertionError,match="Selected Python identity mismatch"):
            exec(code,namespace)
