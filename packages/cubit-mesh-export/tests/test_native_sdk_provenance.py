"""Configured SDK identity must replace stale manifest metadata or fail closed."""
import importlib.util
import json
import shutil
import subprocess
from pathlib import Path

import pytest


ROOT = Path(__file__).resolve().parents[3]
spec = importlib.util.spec_from_file_location(
    "sdk_provenance_test",
    ROOT / "packages/cubit-mesh-export/src/cubit_mesh_export/_native_provenance.py")
provenance = importlib.util.module_from_spec(spec)
spec.loader.exec_module(provenance)


@pytest.mark.parametrize("present", [True, False])
def test_native_cmake_honors_explicit_pybind_config_or_fails(tmp_path, present):
    """Guards: caller-selected dependency config; Assumes: CMake, no compiler/SDK."""
    cmake = shutil.which("cmake")
    assert cmake, "CMake required for native dependency discovery contract"
    source = (ROOT / "src/cubit_plugin/CMakeLists.txt").read_text(encoding="utf-8")
    begin = source.index("# Honor the caller's explicit config directory")
    end = source.index("\nif(pybind11_FOUND", begin)
    # Execute the production discovery block, not a second test implementation.
    config_dir = tmp_path / "explicit-config"
    config_dir.mkdir()
    if present:
        (config_dir / "pybind11Config.cmake").write_text(
            'set(EXPLICIT_CONFIG_LOADED TRUE)\n', encoding="utf-8")
    script = tmp_path / "probe.cmake"
    marker = tmp_path / "discovery.txt"
    script.write_text(source[begin:end] +
        f'\nfile(WRITE "{marker.as_posix()}" "${{EXPLICIT_CONFIG_LOADED}}")\n',
        encoding="utf-8")
    result = subprocess.run([cmake, f"-Dpybind11_DIR={config_dir}",
        "-DCMAKE_FIND_USE_PACKAGE_REGISTRY=FALSE", "-DCMAKE_FIND_USE_SYSTEM_PACKAGE_REGISTRY=FALSE",
        "-DCMAKE_FIND_USE_CMAKE_SYSTEM_PATH=FALSE", "-DCMAKE_FIND_USE_SYSTEM_ENVIRONMENT_PATH=FALSE",
        "-P", str(script)], capture_output=True, text=True, timeout=20)
    assert (result.returncode == 0) is present, (result.returncode, result.stdout, result.stderr)
    if present:
        assert marker.read_text() == "TRUE", "Caller config was not executed"
    else:
        assert "pybind11" in result.stderr, result.stderr
        assert not marker.exists(), "Production discovery did not fail closed"


@pytest.fixture
def candidate(tmp_path, monkeypatch):
    package = tmp_path / "package"
    package.mkdir()
    (package / "native_payloads.json").write_text(json.dumps({
        "payloads": {name: {"cubit_version": "stale"}
                     for name in provenance.REQUIRED_PAYLOADS}}))
    for name in provenance.REQUIRED_PAYLOADS:
        (package / name).write_bytes(b"candidate")
    sdk = tmp_path / "arbitrary-sdk-location"
    (sdk / "cmake").mkdir(parents=True)
    for build in ("build-pyd", "build-ccm"):
        (tmp_path / "src/cubit_plugin" / build).mkdir(parents=True)
    monkeypatch.setattr(provenance, "_source_commit", lambda root: "test")
    monkeypatch.setattr(provenance.metadata, "version", lambda name: "6.2.2607")
    return tmp_path, package, sdk


def configure(root, sdk, cache_type="PATH"):
    for build in ("build-pyd", "build-ccm"):
        (root / "src/cubit_plugin" / build / "CMakeCache.txt").write_text(
            f"Cubit_DIR:{cache_type}={sdk / 'cmake'}\n", encoding="utf-8")


@pytest.mark.parametrize("sdk_version", ["2025.12", "2026.8"])
@pytest.mark.parametrize("cache_type", ["PATH", "UNINITIALIZED"])
def test_record_replaces_stale_sdk_from_build_info(candidate, sdk_version, cache_type):
    root, package, sdk = candidate
    configure(root, sdk, cache_type)
    (sdk / "BUILD_INFO").write_text(f"short={sdk_version}\n", encoding="utf-8")
    manifest = provenance.record_manifest(root, package)
    assert {p["cubit_version"] for p in manifest["payloads"].values()} == {sdk_version}


@pytest.mark.parametrize("build", ["build-pyd", "build-ccm"])
@pytest.mark.parametrize("missing", ["cache", "entry", "empty-entry", "build-info", "short", "empty-short"])
def test_missing_sdk_identity_fails_without_rewriting_manifest(candidate, build, missing):
    root, package, sdk = candidate
    configure(root, sdk)
    info = sdk / "BUILD_INFO"
    info.write_text("short=2026.8\n", encoding="utf-8")
    cache = root / "src/cubit_plugin" / build / "CMakeCache.txt"
    if missing == "cache":
        cache.unlink()
    elif missing in ("entry", "empty-entry"):
        cache.write_text("" if missing == "entry" else "Cubit_DIR:PATH=\n")
    elif missing == "build-info":
        info.unlink()
    else:
        info.write_text("other=value\n" if missing == "short" else "short=\nother=value\n")
    before = (package / "native_payloads.json").read_bytes()
    expected = cache if missing in ("cache", "entry", "empty-entry") else info
    with pytest.raises(ValueError) as error:
        provenance.record_manifest(root, package)
    assert str(expected) in str(error.value)
    assert (package / "native_payloads.json").read_bytes() == before
