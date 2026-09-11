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


def test_release_head_uses_safe_git_helper(monkeypatch):
    calls = []

    def fake_git(*args, **kwargs):
        calls.append((args, kwargs))
        return subprocess.CompletedProcess(args, 0, stdout="ABCDEF\n", stderr="")

    monkeypatch.setattr(release_quad, "_git", fake_git)

    assert release_quad._release_head() == "abcdef"
    assert calls == [(('rev-parse', 'HEAD'), {})]


def test_editable_release_roots_can_target_one_clean_nas_worktree(monkeypatch):
    monkeypatch.setenv(
        release_quad.EDITABLE_REPO_LAB_ENV,
        "S:/Radia/release-quad/Radia-v4.95.46/",
    )
    monkeypatch.setenv(
        release_quad.EDITABLE_REPO_100_ENV,
        "W:\\00_CAE\\Radia\\release-quad\\Radia-v4.95.46\\",
    )

    assert release_quad._lab_editable_packages()[:3] == [
        ("radia", "S:/Radia/release-quad/Radia-v4.95.46"),
        (
            "cubit-mesh-export",
            "S:/Radia/release-quad/Radia-v4.95.46/packages/cubit-mesh-export",
        ),
        (
            "radia-mcp",
            "S:/Radia/release-quad/Radia-v4.95.46/packages/radia-mcp",
        ),
    ]
    assert release_quad._remote_100_editable_packages()[0] == (
        "radia",
        r"W:\00_CAE\Radia\release-quad\Radia-v4.95.46",
    )
    assert release_quad._canonical_lab_editable_packages()[:3] == [
        ("radia", "S:/Radia/01_GitHub"),
        ("cubit-mesh-export", "S:/Radia/01_GitHub/packages/cubit-mesh-export"),
        ("radia-mcp", "S:/Radia/01_GitHub/packages/radia-mcp"),
    ]
    assert release_quad._canonical_remote_100_editable_packages()[0] == (
        "radia",
        r"W:\00_CAE\Radia\01_GitHub",
    )


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


def test_lab_deploy_stops_before_killing_processes_on_source_mismatch(monkeypatch):
    monkeypatch.setattr(release_quad, "_release_head", lambda: "a" * 40)
    monkeypatch.setattr(
        release_quad, "_verify_local_release_source", lambda _repo, _sha: 4
    )

    def unexpected_kill():
        raise AssertionError("processes must not be killed for an invalid source")

    monkeypatch.setattr(release_quad, "_kill_cubit_local", unexpected_kill)
    monkeypatch.setattr(release_quad, "_kill_mcp_local", unexpected_kill)

    assert release_quad._deploy_lab() == 4


def test_remote_deploy_checks_exact_source_before_install(monkeypatch):
    expected_sha = "a" * 40
    captured = {}
    monkeypatch.setattr(release_quad, "_release_head", lambda: expected_sha)

    def capture_run(command, **_kwargs):
        captured["command"] = command

    monkeypatch.setattr(release_quad, "run", capture_run)

    assert release_quad._deploy_editable_remote(
        "release-host", "release host", r"W:\Radia\release-source"
    ) == 0
    script = base64.b64decode(captured["command"][-1]).decode("utf-16le")
    assert expected_sha in script
    assert 'safe.directory=W:/Radia/release-source' in script
    assert "status --porcelain --untracked-files=no" in script
    assert script.index("rev-parse HEAD") < script.index("pip uninstall")
    assert script.index("pip uninstall") < script.index("pip install --no-deps")


def test_remote_restore_forces_canonical_uninstall_then_editable_install(monkeypatch):
    captured = {}

    def capture_run(command, **_kwargs):
        captured["command"] = command
        return subprocess.CompletedProcess(command, 0)

    monkeypatch.setattr(release_quad, "run", capture_run)

    assert release_quad._restore_100_canonical_editable() == 0
    script = base64.b64decode(captured["command"][-1]).decode("utf-16le")
    assert "release-quad" not in script
    assert r"W:\00_CAE\Radia\01_GitHub\packages\radia-mcp" in script
    assert "pip uninstall -y radia cubit-mesh-export radia-mcp" in script
    assert script.index("pip uninstall") < script.index("pip install")
    assert "mcp-server-grant-writing --selftest" in script


def test_done_keeps_exact_verified_editables_after_all_gates(monkeypatch):
    calls = []
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

    def unexpected_restore(_args):
        raise AssertionError("done must not mutate a verified editable tier")

    monkeypatch.setattr(release_quad, "cmd_restore_editable", unexpected_restore)

    args = type("Args", (), {"simulink_package": None})()
    assert release_quad.cmd_done(args) == 0
    assert calls == [
        "preflight", "source", "tag", "lab", "100", "phase9", "guard", "main"
    ]


def test_done_stops_before_machine_checks_when_active_source_is_stale(monkeypatch):
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
