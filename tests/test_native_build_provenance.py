"""Execute the native-build provenance helper against temporary Git repositories."""

from __future__ import annotations

import json
import os
import subprocess
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "tools" / "native_build_provenance.ps1"
BUILD_SCRIPT = ROOT / "Build.ps1"


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


def test_initial_dirty_source_completes_without_manifest(tmp_path):
    repo, binary = _repo(tmp_path)
    (repo / "source.txt").write_text("already dirty\n", encoding="utf-8")
    manifest_path = Path(f"{binary}.build.json")
    command = (
        f". {_quote(SCRIPT)}; "
        f"$start=Get-NativeBuildSourceIdentity -RepoRoot {_quote(repo)}; "
        f"[IO.File]::WriteAllText({_quote(manifest_path)}, 'stale'); "
        f"Write-NativeBuildProvenance -BinaryPath {_quote(binary)} "
        f"-RepoRoot {_quote(repo)} -StartIdentity $start; "
        f"if (Test-Path -LiteralPath {_quote(manifest_path)}) {{ exit 8 }}"
    )
    result = _pwsh(command)
    assert result.returncode == 0, (result.stdout, result.stderr)
    assert "provenance was not issued" in result.stdout
    assert not manifest_path.exists()


def test_dirty_source_change_during_build_is_still_rejected(tmp_path):
    repo, binary = _repo(tmp_path)
    (repo / "source.txt").write_text("dirty at start\n", encoding="utf-8")
    mutation = f"[IO.File]::WriteAllText({_quote(repo / 'source.txt')}, 'dirtier')"
    _assert_change_rejected(repo, binary, mutation)


def test_git_failure_rejects_and_removes_stale_manifest(tmp_path):
    repo, binary = _repo(tmp_path)
    mutation = f"Rename-Item -LiteralPath {_quote(repo / '.git')} -NewName '.git-away'"
    _assert_change_rejected(repo, binary, mutation)


def test_build_entrypoint_declares_optional_strict_provenance_mode():
    source = BUILD_SCRIPT.read_text(encoding="utf-8-sig")
    assert "[switch]$RequireNativeProvenance" in source
    assert "$RequireNativeProvenance -and $NativeBuildSourceIdentity.source_dirty" in source
    assert "Write-NativeBuildProvenance" in source


@pytest.mark.skipif(
    os.environ.get("RADIA_RUN_NATIVE_BUILD_CONTRACT") != "1",
    reason="set RADIA_RUN_NATIVE_BUILD_CONTRACT=1 to exercise real native builds",
)
def test_real_build_clean_certificate_and_dirty_noncertificate():
    # Keep the mapped-drive spelling on Windows. Path.resolve() canonicalizes
    # this checkout to UNC, while Build.ps1's native toolchain requires a drive.
    build_root = Path.cwd()
    assert (build_root / "Build.ps1").samefile(BUILD_SCRIPT)
    build_script = build_root / "Build.ps1"
    status = _git(build_root, "status", "--porcelain=v1", "--untracked-files=all")
    assert not status, f"real build contract requires a clean checkout:\n{status}"
    binary = build_root / "src" / "radia" / "_radia_pybind.pyd"
    manifest = Path(f"{binary}.build.json")
    sentinel = build_root / ".native-build-contract-dirty"

    def run_build(*extra: str) -> subprocess.CompletedProcess[str]:
        environment = os.environ.copy()
        environment["PATH"] = os.pathsep.join(
            [str(Path(sys.executable).parent), environment["PATH"]]
        )
        return subprocess.run(
            [
                "pwsh",
                "-NoLogo",
                "-NoProfile",
                "-ExecutionPolicy",
                "Bypass",
                "-File",
                str(build_script),
                "-RadiaOnly",
                *extra,
            ],
            cwd=build_root,
            text=True,
            capture_output=True,
            check=False,
            env=environment,
            timeout=300,
        )

    clean = run_build("-RequireNativeProvenance")
    assert clean.returncode == 0, (clean.stdout, clean.stderr)
    certificate = json.loads(manifest.read_text(encoding="utf-8"))
    assert certificate["source_dirty"] is False
    assert certificate["source_commit"] == _git(build_root, "rev-parse", "HEAD")

    try:
        sentinel.write_text("intentional dirty-build contract fixture\n", encoding="utf-8")
        dirty = run_build()
        assert dirty.returncode == 0, (dirty.stdout, dirty.stderr)
        assert "provenance was not issued" in dirty.stdout
        assert not manifest.exists()

        strict_dirty = run_build("-RequireNativeProvenance")
        assert strict_dirty.returncode != 0
        assert "source checkout is dirty" in (strict_dirty.stdout + strict_dirty.stderr)
        assert not manifest.exists()
    finally:
        sentinel.unlink(missing_ok=True)
