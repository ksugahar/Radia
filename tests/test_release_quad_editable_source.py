from __future__ import annotations

import base64
import importlib.util
import json
import shutil
import subprocess
import time
from pathlib import Path


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


def test_local_release_source_requires_exact_sha_and_tracked_clean(tmp_path):
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

    assert release_quad._verify_local_release_source(str(repo), head) == 0
    assert release_quad._verify_local_release_source(str(repo), "0" * 40) == 4

    (repo / "tracked.txt").write_text("parallel WIP\n", encoding="ascii")
    assert release_quad._verify_local_release_source(str(repo), head) == 4


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


def test_done_restores_canonical_editables_only_after_all_gates(monkeypatch):
    calls = []
    monkeypatch.setattr(release_quad, "cmd_preflight", lambda _args: calls.append("preflight") or 0)
    monkeypatch.setattr(release_quad, "_verify_lab_editable", lambda *_args: calls.append("lab") or 0)
    monkeypatch.setattr(release_quad, "_verify_100_editable", lambda *_args: calls.append("100") or 0)
    monkeypatch.setattr(release_quad, "cmd_phase9", lambda _args: calls.append("phase9") or 0)
    monkeypatch.setattr(
        release_quad,
        "_run_retired_standalone_pyside_guard",
        lambda: calls.append("guard") or 0,
    )
    monkeypatch.setattr(release_quad, "_check_main_synced", lambda **_kwargs: calls.append("main") or 0)
    monkeypatch.setattr(release_quad, "cmd_restore_editable", lambda _args: calls.append("restore") or 0)

    args = type("Args", (), {"simulink_package": None})()
    assert release_quad.cmd_done(args) == 0
    assert calls == ["preflight", "lab", "100", "phase9", "guard", "main", "restore"]


def test_ci_output_selection_requires_fresh_matching_commit(tmp_path, monkeypatch):
    monkeypatch.setattr(release_quad, "CI_OUTPUT_ROOT", tmp_path)
    started_at = time.time()
    expected_sha = "a" * 40

    wrong = tmp_path / "radia-ci-output-1-1"
    wrong.mkdir()
    (wrong / release_quad.CI_CONTEXT_NAME).write_text(
        json.dumps({
            "schema": "radia.ci-output-context.v1",
            "sha": "b" * 40,
        }),
        encoding="utf-8",
    )

    matching = tmp_path / "radia-ci-output-2-1"
    matching.mkdir()
    (matching / release_quad.CI_CONTEXT_NAME).write_text(
        json.dumps({
            "schema": "radia.ci-output-context.v1",
            "sha": expected_sha,
        }),
        encoding="utf-8",
    )

    assert release_quad._find_ci_output_dir(expected_sha, started_at) == matching
    assert release_quad._find_ci_output_dir("c" * 40, started_at) is None
