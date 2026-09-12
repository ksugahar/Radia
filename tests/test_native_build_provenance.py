"""Execute the native-build provenance helper against temporary Git repositories."""

from __future__ import annotations

import json
import subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "tools" / "native_build_provenance.ps1"


def _git(repo: Path, *arguments: str) -> str:
    return subprocess.check_output(
        ["git", "-C", str(repo), *arguments], text=True
    ).strip()


def _repo(tmp_path: Path) -> tuple[Path, Path]:
    repo = tmp_path / "repo"
    repo.mkdir()
    _git(repo, "init")
    _git(repo, "config", "user.email", "test@example.invalid")
    _git(repo, "config", "user.name", "Provenance Test")
    (repo / "source.txt").write_text("v1\n", encoding="utf-8")
    (repo / ".gitignore").write_text("*.pyd\n*.pyd.build.json\n", encoding="utf-8")
    _git(repo, "add", "source.txt", ".gitignore")
    _git(repo, "commit", "-m", "initial")
    binary = repo / "native.pyd"
    binary.write_bytes(b"native")
    return repo, binary


def _quote(path: Path) -> str:
    return "'" + str(path).replace("'", "''") + "'"


def _pwsh(command: str) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        [
            "pwsh",
            "-NoLogo",
            "-NoProfile",
            "-ExecutionPolicy",
            "Bypass",
            "-Command",
            command,
        ],
        text=True,
        capture_output=True,
        check=False,
    )


def test_dot_source_does_not_change_caller_strict_mode_or_emit_values():
    command = (
        "$ErrorActionPreference='Stop'; "
        f"$loaded=@(. {_quote(SCRIPT)}); "
        "if ($loaded.Count -ne 0) { exit 7 }; "
        "$callerCanReadUndefined=$undefinedAfterDotSource; "
        "if ($null -ne $callerCanReadUndefined) { exit 8 }"
    )
    result = _pwsh(command)
    assert result.returncode == 0, (result.stdout, result.stderr)


def test_manifest_is_emitted_only_when_start_and_end_identity_match(tmp_path):
    repo, binary = _repo(tmp_path)
    manifest_path = Path(f"{binary}.build.json")
    command = (
        f". {_quote(SCRIPT)}; "
        f"$start=Get-NativeBuildSourceIdentity -RepoRoot {_quote(repo)}; "
        f"Write-NativeBuildProvenance -BinaryPath {_quote(binary)} "
        f"-RepoRoot {_quote(repo)} -StartIdentity $start; "
        f"Get-Content -Raw -LiteralPath {_quote(manifest_path)}"
    )
    result = _pwsh(command)
    assert result.returncode == 0, result.stderr
    manifest = json.loads(result.stdout[result.stdout.index("{") :])
    assert manifest["schema"] == "radia.native-build-provenance.v1"
    assert manifest["source_commit"] == _git(repo, "rev-parse", "HEAD")
    assert manifest["source_dirty"] is False
    assert len(manifest["source_change_fingerprint_sha256"]) == 64


def _assert_change_rejected(repo: Path, binary: Path, mutation: str) -> None:
    manifest_path = Path(f"{binary}.build.json")
    command = (
        f". {_quote(SCRIPT)}; "
        f"$start=Get-NativeBuildSourceIdentity -RepoRoot {_quote(repo)}; "
        f"[IO.File]::WriteAllText({_quote(manifest_path)}, 'stale'); "
        f"{mutation}; "
        "try { "
        f"Write-NativeBuildProvenance -BinaryPath {_quote(binary)} "
        f"-RepoRoot {_quote(repo)} -StartIdentity $start; exit 9 "
        "} catch { "
        f"if (Test-Path -LiteralPath {_quote(manifest_path)}) {{ exit 8 }}; exit 0 "
        "}"
    )
    result = _pwsh(command)
    assert result.returncode == 0, (result.stdout, result.stderr)
    assert not manifest_path.exists()


def test_commit_change_during_build_rejects_and_removes_stale_manifest(tmp_path):
    repo, binary = _repo(tmp_path)
    mutation = (
        f"[IO.File]::WriteAllText({_quote(repo / 'source.txt')}, 'v2'); "
        f"git -C {_quote(repo)} add source.txt; "
        f"git -C {_quote(repo)} commit -m changed | Out-Null"
    )
    _assert_change_rejected(repo, binary, mutation)


def test_dirty_state_change_during_build_rejects_and_removes_stale_manifest(tmp_path):
    repo, binary = _repo(tmp_path)
    mutation = f"[IO.File]::WriteAllText({_quote(repo / 'source.txt')}, 'dirty')"
    _assert_change_rejected(repo, binary, mutation)


def test_initial_dirty_source_never_emits_manifest(tmp_path):
    repo, binary = _repo(tmp_path)
    (repo / "source.txt").write_text("already dirty\n", encoding="utf-8")
    manifest_path = Path(f"{binary}.build.json")
    command = (
        f". {_quote(SCRIPT)}; "
        f"$start=Get-NativeBuildSourceIdentity -RepoRoot {_quote(repo)}; "
        f"[IO.File]::WriteAllText({_quote(manifest_path)}, 'stale'); "
        "try { "
        f"Write-NativeBuildProvenance -BinaryPath {_quote(binary)} "
        f"-RepoRoot {_quote(repo)} -StartIdentity $start; exit 9 "
        "} catch { "
        f"if (Test-Path -LiteralPath {_quote(manifest_path)}) {{ exit 8 }}; exit 0 "
        "}"
    )
    result = _pwsh(command)
    assert result.returncode == 0, (result.stdout, result.stderr)
    assert not manifest_path.exists()


def test_git_failure_rejects_and_removes_stale_manifest(tmp_path):
    repo, binary = _repo(tmp_path)
    mutation = f"Rename-Item -LiteralPath {_quote(repo / '.git')} -NewName '.git-away'"
    _assert_change_rejected(repo, binary, mutation)
