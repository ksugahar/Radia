from __future__ import annotations

import base64
import importlib.util
import subprocess
from pathlib import Path


_TOOL = Path(__file__).resolve().parents[1] / "tools" / "release_qud.py"
_SPEC = importlib.util.spec_from_file_location("radia_release_qud_tool", _TOOL)
assert _SPEC is not None and _SPEC.loader is not None
release_qud = importlib.util.module_from_spec(_SPEC)
_SPEC.loader.exec_module(release_qud)


def _git(repo, *args):
    return subprocess.run(
        ["git", "-C", str(repo), *args],
        check=True,
        capture_output=True,
        text=True,
    ).stdout.strip()


def test_editable_release_roots_can_target_one_clean_nas_worktree(monkeypatch):
    monkeypatch.setenv(
        release_qud.EDITABLE_REPO_LAB_ENV,
        "S:/Radia/release-qud/Radia-v4.95.46/",
    )
    monkeypatch.setenv(
        release_qud.EDITABLE_REPO_100_ENV,
        "W:\\00_CAE\\Radia\\release-qud\\Radia-v4.95.46\\",
    )

    assert release_qud._lab_editable_packages()[:3] == [
        ("radia", "S:/Radia/release-qud/Radia-v4.95.46"),
        (
            "cubit-mesh-export",
            "S:/Radia/release-qud/Radia-v4.95.46/packages/cubit-mesh-export",
        ),
        (
            "radia-mcp",
            "S:/Radia/release-qud/Radia-v4.95.46/packages/radia-mcp",
        ),
    ]
    assert release_qud._remote_100_editable_packages()[0] == (
        "radia",
        r"W:\00_CAE\Radia\release-qud\Radia-v4.95.46",
    )


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

    assert release_qud._verify_local_release_source(str(repo), head) == 0
    assert release_qud._verify_local_release_source(str(repo), "0" * 40) == 4

    (repo / "tracked.txt").write_text("parallel WIP\n", encoding="ascii")
    assert release_qud._verify_local_release_source(str(repo), head) == 4


def test_lab_deploy_stops_before_killing_processes_on_source_mismatch(monkeypatch):
    monkeypatch.setattr(release_qud, "_release_head", lambda: "a" * 40)
    monkeypatch.setattr(
        release_qud, "_verify_local_release_source", lambda _repo, _sha: 4
    )

    def unexpected_kill():
        raise AssertionError("processes must not be killed for an invalid source")

    monkeypatch.setattr(release_qud, "_kill_cubit_local", unexpected_kill)
    monkeypatch.setattr(release_qud, "_kill_mcp_local", unexpected_kill)

    assert release_qud._deploy_lab() == 4


def test_remote_deploy_checks_exact_source_before_install(monkeypatch):
    expected_sha = "a" * 40
    captured = {}
    monkeypatch.setattr(release_qud, "_release_head", lambda: expected_sha)

    def capture_run(command, **_kwargs):
        captured["command"] = command

    monkeypatch.setattr(release_qud, "run", capture_run)

    assert release_qud._deploy_editable_remote(
        "release-host", "release host", r"W:\Radia\release-source"
    ) == 0
    script = base64.b64decode(captured["command"][-1]).decode("utf-16le")
    assert expected_sha in script
    assert 'safe.directory=W:/Radia/release-source' in script
    assert "status --porcelain --untracked-files=no" in script
    assert script.index("rev-parse HEAD") < script.index("pip install -e")
