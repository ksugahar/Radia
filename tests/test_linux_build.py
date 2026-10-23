"""Exercise HPC build admission checks and failure propagation without compiling."""

import importlib.util
import json
from pathlib import Path
import subprocess

import pytest

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location("linux_builder", ROOT / "tools/build_linux.py")
builder = importlib.util.module_from_spec(spec)
spec.loader.exec_module(builder)


def test_versions_use_package_metadata_and_reject_mismatch(monkeypatch):
    pins = {"ngsolve": "6.2.2608", "netgen-mesher": "6.2.2608"}
    monkeypatch.setattr(builder.importlib.metadata, "version", pins.__getitem__)
    assert builder.check_versions() == pins
    pins["netgen-mesher"] = "6.2.2607"
    with pytest.raises(RuntimeError, match="netgen-mesher: installed 6.2.2607, required 6.2.2608"):
        builder.check_versions()


@pytest.mark.parametrize("failure", ["missing-compiler", "configure", "build", "import", None])
def test_pipeline_stops_on_failure_and_records_artifacts(tmp_path, monkeypatch, failure):
    source = tmp_path / "checkout"
    artifacts = source / "src/radia"
    artifacts.mkdir(parents=True)
    artifact = artifacts / "_radia_pybind.so"
    artifact.write_bytes(b"fixture native module")
    mkl = tmp_path / "mkl"
    (mkl / "include").mkdir(parents=True)
    (mkl / "lib").mkdir()
    (mkl / "include/mkl.h").touch()
    (mkl / "lib/libmkl_rt.so").touch()
    monkeypatch.setenv("MKLROOT", str(mkl))
    monkeypatch.delenv("CC", raising=False)
    monkeypatch.delenv("CXX", raising=False)
    monkeypatch.setattr(builder, "ROOT", source)
    monkeypatch.setattr(builder.sys, "platform", "linux")
    monkeypatch.setattr(builder.sys, "version_info", (3, 12))
    monkeypatch.setattr(builder.os, "sched_getaffinity", lambda _: {20}, raising=False)
    monkeypatch.setattr(builder, "check_versions", lambda: {"ngsolve": "6.2.2608"})
    monkeypatch.setattr(builder.shutil, "which", lambda name:
                        None if name == "icpx" and failure == "missing-compiler" else name)
    commands = []

    def run(command, **kwargs):
        command = list(map(str, command))
        commands.append(command)
        stage = ("configure" if "-S" in command else "build" if "--build" in command
                 else "import" if "-c" in command else "other")
        if stage == failure:
            raise subprocess.CalledProcessError(7, command)
        output = b"diff" if kwargs.get("capture_output") and not kwargs.get("text") else ""
        if "--version" in command:
            output = "Intel(R) oneAPI DPC++/C++ Compiler"
        if "rev-parse" in command:
            output = "a" * 40
        return subprocess.CompletedProcess(command, 0, output, "")

    monkeypatch.setattr(builder, "run", run)
    build = tmp_path / "build"
    result = builder.main(["--build-dir", str(build), "--targets", "_radia_pybind"])
    report = json.loads((build / "linux-build.json").read_text())
    assert result == (1 if failure else 0), report
    assert report["passed"] == (failure is None), report
    configured = [c for c in commands if "-S" in c]
    built = [c for c in commands if "--build" in c]
    imported = [c for c in commands if "-c" in c]
    if failure == "missing-compiler":
        assert not configured and not built and not imported, commands
    elif failure == "configure":
        assert configured and not built and not imported, commands
    elif failure == "build":
        assert built and not imported, commands
    else:
        assert imported, commands
    if configured:
        assert "-DCMAKE_CXX_COMPILER=icpx" in configured[0], configured
    if built:
        assert built[0][-2:] == ["--parallel", "1"], built
    if failure is None:
        import hashlib
        assert report["artifacts"]["_radia_pybind"]["sha256"] == hashlib.sha256(artifact.read_bytes()).hexdigest()
        manifest = json.loads(Path(str(artifact) + ".build.json").read_text())
        assert manifest["source_commit"] == "a" * 40
        assert manifest["binary_sha256"] == report["artifacts"]["_radia_pybind"]["sha256"]
    else:
        assert report["artifacts"] == {}, report


def test_dirty_source_cannot_receive_acceptance_manifest(tmp_path):
    binary = tmp_path / "native.so"
    binary.write_bytes(b"native")
    assert builder.write_provenance(binary, {"source_dirty": True}) is False
    assert not Path(str(binary) + ".build.json").exists()


def test_source_identity_matches_windows_provenance_helper(tmp_path):
    """Use the existing Windows implementation as an independent contract oracle."""
    import shutil
    if shutil.which("pwsh") is None:
        pytest.skip("PowerShell is needed for the Windows/Linux contract comparison")
    repo = tmp_path / "repo"
    repo.mkdir()
    def git(*args):
        subprocess.run(["git", "-C", str(repo), *args], check=True, capture_output=True)
    git("init")
    git("config", "user.name", "Build Contract Test")
    git("config", "user.email", "test@example.invalid")
    (repo / "source.txt").write_text("original\n", encoding="utf-8")
    git("add", "source.txt")
    git("commit", "-m", "fixture")
    helper = str(ROOT / "tools/native_build_provenance.ps1").replace("'", "''")
    quoted_repo = str(repo).replace("'", "''")
    for content in ("original\n", "changed\n"):
        (repo / "source.txt").write_text(content, encoding="utf-8")
        windows = subprocess.run(
            ["pwsh", "-NoProfile", "-Command",
             f". '{helper}'; Get-NativeBuildSourceIdentity -RepoRoot '{quoted_repo}' | ConvertTo-Json -Compress"],
            check=True, capture_output=True, text=True, encoding="utf-8")
        assert builder.source_identity(repo) == json.loads(windows.stdout)
