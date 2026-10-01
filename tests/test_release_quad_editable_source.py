from __future__ import annotations

import base64
import importlib.util
import shutil
import subprocess
from pathlib import Path

import pytest


_TOOL = Path(__file__).resolve().parents[1] / "tools" / "release_quad.py"
_SPEC = importlib.util.spec_from_file_location("radia_release_quad_tool", _TOOL)
assert _SPEC is not None and _SPEC.loader is not None
release_quad = importlib.util.module_from_spec(_SPEC)
_SPEC.loader.exec_module(release_quad)
_GIT = shutil.which("git")
assert _GIT is not None, "Git is required by the release-quad contract tests"


def _git(repo, *args):
    return subprocess.run(
        [_GIT, "-C", str(repo), *args],
        check=True,
        capture_output=True,
        text=True,
    ).stdout.strip()


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
    calls = []

    def fake_git(*args, **kwargs):
        calls.append((args, kwargs))
        return subprocess.CompletedProcess(
            args, 0, stdout="a" * 40 + "\n", stderr="")

    monkeypatch.setattr(release_quad, "_git", fake_git)

    assert release_quad._release_head() == "a" * 40
    assert calls == [(('rev-list', '-n', '1', 'v5.1.0'), {'check': False})]


def test_phase8_refuses_editable_targets_from_100_before_install(monkeypatch):
    from argparse import Namespace

    monkeypatch.setattr(release_quad, "cmd_preflight", lambda _args: 0)
    monkeypatch.setattr(release_quad.platform, "node", lambda: "INTEL11")
    monkeypatch.setattr(
        release_quad, "_deploy_lab",
        lambda: pytest.fail("LAB install must not run from 100号機"),
    )
    monkeypatch.setattr(
        release_quad, "_deploy_100",
        lambda: pytest.fail("editable deployment must use the LAB controller"),
    )

    assert release_quad.cmd_phase8(Namespace(target="lab,100")) == 2


def test_phase8_routes_editable_targets_from_lab(monkeypatch):
    from argparse import Namespace

    calls = []
    monkeypatch.setattr(release_quad, "cmd_preflight", lambda _args: 0)
    monkeypatch.setattr(release_quad.platform, "node", lambda: "LAB")
    monkeypatch.setattr(release_quad, "_deploy_lab", lambda: calls.append("lab") or 0)
    monkeypatch.setattr(release_quad, "_deploy_100", lambda: calls.append("100") or 0)

    assert release_quad.cmd_phase8(Namespace(target="lab,100")) == 0
    assert calls == ["lab", "100"]


def test_phase8_all_refuses_non_lab_controller(monkeypatch):
    from argparse import Namespace

    monkeypatch.setattr(release_quad, "cmd_preflight", lambda _args: 0)
    monkeypatch.setattr(release_quad.platform, "node", lambda: "INTEL11")
    monkeypatch.setattr(
        release_quad, "_deploy_lab",
        lambda: pytest.fail("all must fail before its first install"),
    )

    assert release_quad.cmd_phase8(Namespace(target="all")) == 2


def test_default_editable_roots_are_fixed_release_checkouts(monkeypatch):
    monkeypatch.delenv(release_quad.EDITABLE_REPO_LAB_ENV, raising=False)
    monkeypatch.delenv(release_quad.EDITABLE_REPO_100_ENV, raising=False)
    monkeypatch.setattr(release_quad, "_radia_version", lambda: "5.1.0")
    monkeypatch.setattr(release_quad, "_release_commit", lambda: "c" * 40)

    assert release_quad._editable_repo_lab() == (
        "S:/Radia/release-quad/v5.1.0-ccccccccc")
    assert release_quad._editable_repo_100() == (
        r"W:\00_CAE\Radia\release-quad\v5.1.0-ccccccccc")


def test_editable_runtime_gate_requires_exact_native_manifest():
    gate = release_quad.EDITABLE_RELEASE_VERIFY
    assert "release_native_payloads.json" in gate
    assert 'radia.release-native-payloads.v1' in gate
    assert '{"_radia_pybind.pyd", "axifem.pyd", "sparsesolv_ngsolve.pyd"}' in gate
    assert "native payload set differs" in gate
    assert 'name == "cln_core.pyd"' in gate
    assert '"locked-old" in lowered' in gate


def test_deployment_plan_is_solver_only_and_does_not_probe_runtime(
        monkeypatch, capsys):
    import json

    monkeypatch.setattr(release_quad, "run", lambda *a, **k: pytest.fail("not a dry run"))
    monkeypatch.setattr(release_quad, "_pip_show", lambda *a: pytest.fail("runtime inference"))
    assert release_quad.cmd_deployment_plan(None) == 0
    plans = json.loads(capsys.readouterr().out)
    assert [p["solver_source"] for p in plans] == [
        release_quad._editable_repo_lab(), release_quad._editable_repo_100()]
    assert all(p["solver_action"] == "verify-and-install" for p in plans)
    assert all(p["independent_packages"] == "unchanged" for p in plans)
    assert all(not p["verified"] for p in plans)
    assert all("mcp" not in key and "cubit" not in key
               for plan in plans for key in plan)


@pytest.mark.parametrize("drift", [0, 1])
def test_lab_deploy_changes_only_radia(monkeypatch, drift):
    root = release_quad._editable_repo_lab()
    monkeypatch.setattr(release_quad, "_release_head", lambda: "a" * 40)
    monkeypatch.setattr(release_quad, "_verify_local_release_source", lambda *a: 0)
    events = []

    def verify(packages):
        assert packages == [("radia", root)]
        events.append("verify")
        return drift

    monkeypatch.setattr(release_quad, "_verify_lab_editable", verify)
    def run(command, **_kwargs):
        events.append(command)
        return subprocess.CompletedProcess(command, 0)
    monkeypatch.setattr(release_quad, "run", run)
    monkeypatch.setattr(release_quad, "_record_release_intent_lab",
                        lambda repo: events.append(("record", repo)) or 0)

    assert release_quad._deploy_lab() == (4 if drift else 0)
    abi_probe = release_quad._solver_abi_probe_command()
    assert abi_probe in events
    command = next(event for event in events if isinstance(event, list) and "pip" in event)
    assert events.index(abi_probe) < events.index(command)
    joined = " ".join(command)
    assert "pip install" in joined
    assert "uninstall" not in joined
    assert "radia-mcp" not in joined
    assert "cubit-mesh-export" not in joined
    assert "Stop-Process" not in joined
    if drift:
        assert events[-1] == "verify"
        assert not any(isinstance(e, tuple) for e in events)
    else:
        assert events[-2] == "verify"
        assert events[-1] == ("record", root)


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

    assert release_quad._deploy_editable_remote("100", "100", "W:/release") == (4 if drift else 0)
    script = base64.b64decode(calls[0][-1]).decode("utf-16le")
    assert "pip install --no-deps" in script
    assert "pip uninstall" not in script
    assert "radia-mcp" not in script
    assert "cubit-mesh-export" not in script
    assert "Stop-Process" not in script
    assert "status --porcelain --untracked-files=no" in script
    assert "ngsolve.__version__" in script
    assert '& "python" "-c"' in script
    assert '\n"python" "-c"' not in script
    assert "netgen-mesher" in script
    assert script.index("rev-parse HEAD") < script.index("pip install")
    assert script.index("ngsolve.__version__") < script.index("pip install")
    if drift:
        assert calls[-1] == "verify"
        assert not any(isinstance(c, tuple) for c in calls)
    else:
        assert calls[-2] == "verify"
        assert calls[-1] == ("record", "100", "W:/release")


def test_done_checks_only_solver_editable_roots(monkeypatch, tmp_path):
    from argparse import Namespace

    # `done` takes its source from the recorded editable intent and refuses
    # when there is neither a record nor a release override -- no default tree
    # is assumed. Record one in an isolated file so this test measures which
    # roots are checked, not what this machine happens to have recorded.
    intent_file = tmp_path / "editable-intent.json"
    monkeypatch.setenv(release_quad.editable_intent.INTENT_FILE_ENV,
                       str(intent_file))
    monkeypatch.delenv(release_quad.EDITABLE_REPO_LAB_ENV, raising=False)
    monkeypatch.delenv(release_quad.EDITABLE_REPO_100_ENV, raising=False)
    intent_module = release_quad.editable_intent
    data = intent_module.load_intent(intent_file)
    intent_module.set_entry(data, "radia", {
        "source": release_quad._editable_repo_lab(), "commit": None,
        "tracked_clean": None, "recorded_at": "2026-09-21T00:00:00Z",
        "recorded_by": "test", "recorded_via": "test", "reason": "fixture",
        "pushed_refs": None, "previous": None})
    intent_module.save_intent(data, intent_file)

    calls = []
    monkeypatch.setattr(release_quad, "cmd_preflight", lambda a: 0)
    # the retired-override gate reaches mdx1/mdx2/hibino over SSH; stub it
    # like every other gate so this test still measures the editable roots
    monkeypatch.setattr(release_quad, "cmd_temp_shadows", lambda a: 0)
    monkeypatch.setattr(release_quad, "_release_head", lambda: "c" * 40)
    monkeypatch.setattr(release_quad, "_verify_local_release_source", lambda root, sha: calls.append((root, sha)) or 0)
    monkeypatch.setattr(release_quad, "_verify_head_release_tag", lambda: calls.append("tag") or 0)
    monkeypatch.setattr(release_quad, "_verify_lab_editable", lambda: calls.append(dict(release_quad._lab_editable_packages())) or 0)
    monkeypatch.setattr(release_quad, "_verify_100_editable", lambda: calls.append(dict(release_quad._remote_100_editable_packages())) or 0)
    monkeypatch.setattr(release_quad, "cmd_phase9", lambda a: 4)
    assert release_quad.cmd_done(Namespace(simulink_package=None)) == 4
    assert calls[0] == (release_quad._editable_repo_lab(), "c" * 40)
    assert calls[1] == "tag"
    assert calls[2] == {"radia": release_quad._editable_repo_lab()}
    assert calls[3] == {"radia": release_quad._editable_repo_100()}


def test_done_requires_simulink_candidate_for_5_1_and_newer(monkeypatch, tmp_path):
    from argparse import Namespace

    monkeypatch.setattr(release_quad, "cmd_preflight", lambda _a: 0)
    monkeypatch.setattr(release_quad, "_done_active_lab_source", lambda: str(tmp_path))
    monkeypatch.setattr(release_quad, "_verify_local_release_source", lambda *_a: 0)
    monkeypatch.setattr(release_quad, "cmd_temp_shadows", lambda _a: 0)
    monkeypatch.setattr(release_quad, "_verify_head_release_tag", lambda: 0)
    monkeypatch.setattr(release_quad, "_verify_lab_editable", lambda: 0)
    monkeypatch.setattr(release_quad, "_verify_100_editable", lambda: 0)
    monkeypatch.setattr(release_quad, "cmd_phase9", lambda _a: 0)
    monkeypatch.setattr(release_quad, "_run_retired_standalone_pyside_guard", lambda: 0)
    monkeypatch.setattr(release_quad, "_check_main_synced", lambda **_k: 0)
    monkeypatch.setattr(release_quad, "_radia_version", lambda: "5.1.0")

    assert release_quad.cmd_done(Namespace(simulink_package=None)) == 4


@pytest.mark.parametrize("local_rc,remote_rc,expected", [(0, 0, 0), (1, 0, 4), (0, 1, 4)])
def test_done_pip_check_requires_both_editable_hosts(monkeypatch, local_rc, remote_rc, expected):
    calls = []

    def fake_run(command, **_kwargs):
        calls.append(command)
        rc = remote_rc if command[0] == "ssh" else local_rc
        return subprocess.CompletedProcess(command, rc)

    monkeypatch.setattr(release_quad, "run", fake_run)
    assert release_quad._verify_final_pip_checks() == expected
    assert calls[0][-2:] == ["pip", "check"]
    remote = base64.b64decode(calls[1][-1]).decode("utf-16le")
    assert remote == "python -m pip check"


def test_restore_editable_is_a_tombstone_that_restores_nothing(monkeypatch, capsys):
    """The command stays only to refuse and name its replacement."""
    events = []
    monkeypatch.setattr(release_quad, "run",
                        lambda *a, **k: events.append(a) or None)
    assert release_quad.cmd_restore_editable(None) == 2
    out = capsys.readouterr().out
    assert "repoint" in out
    assert not events


def test_editable_release_roots_can_target_one_clean_nas_worktree(monkeypatch):
    monkeypatch.setenv(
        release_quad.EDITABLE_REPO_LAB_ENV,
        "S:/Radia/release-quad/Radia-v4.95.46/",
    )
    monkeypatch.setenv(
        release_quad.EDITABLE_REPO_100_ENV,
        "W:\\00_CAE\\Radia\\release-quad\\Radia-v4.95.46\\",
    )

    assert release_quad._lab_editable_packages() == [
        ("radia", "S:/Radia/release-quad/Radia-v4.95.46")]
    assert release_quad._remote_100_editable_packages() == [
        ("radia", r"W:\00_CAE\Radia\release-quad\Radia-v4.95.46")]


def test_unc_normalization_covers_canonical_and_release_worktrees():
    release_unc = (
        r"\\192.168.121.100\work\00_CAE\Radia\release-qud"
        r"\radia-4.95.75\src\radia\__init__.py"
    )
    legacy_root_unc = (
        r"\\192.168.11.100\work\00_CAE\Radia\01_GitHub"
        r"\packages\radia-mcp"
    )

    assert release_quad._norm_path(release_unc) == (
        "s:/radia/release-qud/radia-4.95.75/src/radia/__init__.py"
    )
    assert release_quad._norm_path(legacy_root_unc) == (
        "s:/radia/01_github/packages/radia-mcp"
    )
    assert release_quad._norm_path(
        r"W:\00_CAE\Radia\01_GitHub\src\radia\__init__.py"
    ) == "s:/radia/01_github/src/radia/__init__.py"
    assert "//192.168.121.100/work/00_cae/radia/" in (
        release_quad.REMOTE_EDITABLE_VERIFY
    )


def test_editable_probe_resolves_git_inside_the_target_process():
    probe = release_quad.CROSS_MACHINE_PROBE_LAB

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


def test_lab_deploy_stops_before_install_on_source_mismatch(monkeypatch):
    monkeypatch.setattr(release_quad, "_release_head", lambda: "a" * 40)
    monkeypatch.setattr(
        release_quad, "_verify_local_release_source", lambda _repo, _sha: 4
    )

    monkeypatch.setattr(
        release_quad, "run",
        lambda *_a, **_k: pytest.fail("install must not run for invalid source"))

    assert release_quad._deploy_lab() == 4


def test_lab_deploy_preserves_install_when_binary_preflight_fails(monkeypatch):
    monkeypatch.setattr(release_quad, "_release_head", lambda: "a" * 40)
    monkeypatch.setattr(release_quad, "_verify_local_release_source", lambda *_: 0)
    calls = []
    def blocked(command, **_kwargs):
        calls.append(command)
        rc = 0 if command == release_quad._solver_abi_probe_command() else 1
        return subprocess.CompletedProcess(command, rc)
    monkeypatch.setattr(release_quad, "run", blocked)
    assert release_quad._deploy_lab() == 3
    assert len(calls) == 2
    assert calls[1][-1] == release_quad.SOLVER_INSTALL_GUARD
    assert not any("pip" in command for command in calls)


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
        "release-host", "release host", r"W:\Radia\release-source"
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


def test_done_keeps_exact_verified_editables_after_all_gates(monkeypatch, tmp_path):
    monkeypatch.setenv(release_quad.EDITABLE_REPO_LAB_ENV, str(tmp_path))
    calls = []
    monkeypatch.setattr(release_quad, "cmd_temp_shadows", lambda _args: calls.append("shadows") or 0)
    monkeypatch.setattr(release_quad, "cmd_preflight", lambda _args: calls.append("preflight") or 0)
    monkeypatch.setattr(release_quad, "_release_head", lambda: "a" * 40)
    monkeypatch.setattr(
        release_quad,
        "_verify_local_release_source",
        lambda _repo, _sha: calls.append("source") or 0,
    )
    monkeypatch.setattr(
        release_quad,
        "_verify_head_release_tag",
        lambda: calls.append("tag") or 0,
    )
    monkeypatch.setattr(release_quad, "_verify_lab_editable", lambda *_args: calls.append("lab") or 0)
    monkeypatch.setattr(release_quad, "_verify_100_editable", lambda *_args: calls.append("100") or 0)
    monkeypatch.setattr(release_quad, "cmd_phase9", lambda _args: calls.append("phase9") or 0)
    monkeypatch.setattr(
        release_quad,
        "_run_retired_standalone_pyside_guard",
        lambda: calls.append("guard") or 0,
    )
    monkeypatch.setattr(release_quad, "_check_main_synced", lambda **_kwargs: calls.append("main") or 0)
    monkeypatch.setattr(
        release_quad, "_verify_simulink_candidate_state",
        lambda package: calls.append(("simulink", package)) or 0)
    monkeypatch.setattr(
        release_quad, "_verify_final_pip_checks",
        lambda: calls.append("pip-check") or 0)

    args = type("Args", (), {"simulink_package": "candidate.zip"})()
    assert release_quad.cmd_done(args) == 0
    assert calls == [
        "preflight", "source", "shadows", "tag", "lab", "100", "phase9", "guard", "main",
        ("simulink", "candidate.zip"), "pip-check"
    ]


def test_done_stops_before_machine_checks_when_active_source_is_stale(monkeypatch, tmp_path):
    monkeypatch.setenv(release_quad.EDITABLE_REPO_LAB_ENV, str(tmp_path))
    calls = []
    monkeypatch.setattr(release_quad, "cmd_preflight", lambda _args: calls.append("preflight") or 0)
    monkeypatch.setattr(release_quad, "_release_head", lambda: "a" * 40)
    monkeypatch.setattr(
        release_quad,
        "_verify_local_release_source",
        lambda _repo, _sha: calls.append("source") or 4,
    )

    def unexpected_machine_check(*_args, **_kwargs):
        raise AssertionError("machine checks must not run from a stale source")

    monkeypatch.setattr(release_quad, "_verify_lab_editable", unexpected_machine_check)
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
