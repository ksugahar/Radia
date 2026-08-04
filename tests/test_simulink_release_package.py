from pathlib import Path
import zipfile

import pytest


ROOT = Path(__file__).resolve().parents[1]


def test_package_builder_requires_native_ih_assets():
    module = load_module(
        "package_simulink_release_assets",
        ROOT / "tools" / "package_simulink_release.py",
    )
    assert module.REQUIRED_MEX == ("radia_mex.mexw64",)
    assert set(module.REQUIRED_MATLAB_SFUNCTIONS) == {
        "radia_ih_eddy_sfun.m",
        "radia_ih_thermal_sfun.m",
        "+radia/+simulink/ihEddySFunction.m",
        "+radia/+simulink/ihThermalSFunction.m",
    }
    assert "radia_ih.slx" in module.REQUIRED_MODELS
    assert any(
        name.endswith("validateIHNativeConfig.m")
        for name in module.PACKAGE_FILES
    )
    assert {
        "+radia/+simulink/addIHGeometryUpdateBlock.m",
        "+radia/+simulink/browseIHGeometryFile.m",
        "+radia/+simulink/fileFingerprint.m",
        "+radia/+simulink/normalizeIHGeometryRoles.m",
        "+radia/+simulink/updateIHGeometry.m",
    } <= set(module.PACKAGE_FILES)
    assert not any(name.startswith("radia_ih_") for name in module.REQUIRED_MEX)


def test_package_builder_fails_when_mex_is_missing(tmp_path):
    import importlib.util

    path = ROOT / "tools" / "package_simulink_release.py"
    spec = importlib.util.spec_from_file_location("package_simulink_release", path)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    with pytest.raises(FileNotFoundError, match="Required native release binary"):
        module.build_package(tmp_path, tmp_path / "out")


def test_package_commit_uses_github_sha_without_git(monkeypatch):
    module = load_module(
        "package_simulink_release_commit",
        ROOT / "tools" / "package_simulink_release.py",
    )
    expected = "A" * 40
    monkeypatch.setenv("GITHUB_SHA", expected)
    monkeypatch.setenv("PATH", "")
    assert module.commit() == expected.lower()


def test_package_is_hashed_native_ih_allowlist(tmp_path):
    package_module = load_module(
        "package_simulink_release", ROOT / "tools" / "package_simulink_release.py"
    )
    verify_module = load_module(
        "verify_simulink_release", ROOT / "tools" / "verify_simulink_release.py"
    )
    mex_dir = tmp_path / "mex"
    mex_dir.mkdir()
    for name in (
        *package_module.REQUIRED_MEX,
        *package_module.FULL_RUNTIME_DLLS,
    ):
        (mex_dir / name).write_bytes(fake_x64_pe())
    archive, sums = package_module.build_package(mex_dir, tmp_path / "out")
    manifest = verify_module.verify_archive(archive)
    assert sums.read_text(encoding="ascii").split()[0] == \
        package_module.sha256(archive)
    assert manifest["release_channel"] == "preview"
    assert manifest["schema"] == "radia.simulink.ih-release-manifest.v2"
    assert manifest["backend"] == "matlab-level2+radia-mex-handles"
    assert manifest["required_mex"] == ["matlab/radia_mex.mexw64"]
    assert manifest["standalone_mex_debug_api"] is True
    assert manifest["python_runtime_required_for_native_mex"] is True
    assert manifest["python_per_step"] is False
    assert manifest["operator_assembly"] == "preassembled"
    assert (tmp_path / "out" / "manifest.json").is_file()
    assert len(sums.read_text(encoding="ascii").splitlines()) == 2
    with zipfile.ZipFile(archive) as bundle:
        names = set(bundle.namelist())
    assert "matlab/radia_ih.slx" in names
    assert "matlab/radia_ih_eddy_sfun.m" in names
    assert "matlab/+radia/+simulink/ihEddySFunction.m" in names
    assert "matlab/radia_ih_eddy_sfun.mexw64" not in names
    assert not any("LUT" in name or "makeIHPlant" in name for name in names)
    assert "matlab/radia_simulink_library.slx" not in names


def test_full_library_package_includes_mex_models_and_runtime(tmp_path):
    package_module = load_module(
        "package_simulink_full_release",
        ROOT / "tools" / "package_simulink_release.py",
    )
    verify_module = load_module(
        "verify_simulink_full_release",
        ROOT / "tools" / "verify_simulink_release.py",
    )
    mex_dir = tmp_path / "mex"
    mex_dir.mkdir()
    for name in (
        *package_module.FULL_REQUIRED_MEX,
        *package_module.FULL_RUNTIME_DLLS,
    ):
        (mex_dir / name).write_bytes(fake_x64_pe())
    archive, sums = package_module.build_package(
        mex_dir,
        tmp_path / "out",
        full_library=True,
    )
    manifest = verify_module.verify_archive(archive)
    assert manifest["schema"] == (
        "radia.simulink.library-release-manifest.v2"
    )
    assert manifest["release_channel"] == "production"
    assert manifest["entry_model"] == "matlab/radia_simulink_library.slx"
    assert manifest["python_per_step"] is False
    assert manifest["python_fallback_per_step"] is False
    assert manifest["backend"] == "application-specific"
    assert manifest["ih_backend"] == "matlab-level2+radia-mex-handles"
    assert manifest["required_mex"] == ["matlab/radia_mex.mexw64"]
    assert manifest["required_matlab_products"] == ["MATLAB", "Simulink"]
    assert manifest["feature_toolbox_requirements"] == {
        "adjoint_topology_optimization": ["Optimization Toolbox"],
        "stream_function_topology_optimization": ["Optimization Toolbox"],
    }
    assert manifest["application_batch_backend"] == (
        "python-headless-or-native-as-declared-by-block"
    )
    assert sums.read_text(encoding="ascii").split()[0] == \
        package_module.sha256(archive)
    with zipfile.ZipFile(archive) as bundle:
        names = set(bundle.namelist())
    assert verify_module.FULL_REQUIRED_MEMBERS <= names
    assert "matlab/+radia/+simulink/motorAngleFamilyMexSFunction.m" in names
    assert "matlab/python_api_parity_manifest.json" in names


def test_full_library_enumerates_only_tracked_matlab_files(monkeypatch, tmp_path):
    package_module = load_module(
        "package_simulink_tracked_files",
        ROOT / "tools" / "package_simulink_release.py",
    )
    tracked = tmp_path / "matlab" / "tracked.m"
    tracked.parent.mkdir()
    tracked.write_text("disp('tracked')\n", encoding="utf-8")
    (tracked.parent / "untracked.m").write_text(
        "disp('untracked')\n", encoding="utf-8"
    )
    calls = []

    class Result:
        stdout = "matlab/tracked.m\n"

    def fake_run(command, **kwargs):
        calls.append((command, kwargs))
        return Result()

    monkeypatch.setattr(package_module, "ROOT", tmp_path)
    monkeypatch.setattr(package_module.subprocess, "run", fake_run)

    assert package_module.release_matlab_files() == (Path("matlab/tracked.m"),)
    assert calls[0][0] == [
        "git",
        "-c",
        f"safe.directory={tmp_path.as_posix()}",
        "ls-files",
        "--cached",
        "--",
        "matlab",
    ]


def test_matlab_smoke_decodes_utf8_without_cp932(monkeypatch, tmp_path):
    package_module = load_module(
        "package_simulink_release_encoding",
        ROOT / "tools" / "package_simulink_release.py",
    )
    verify_module = load_module(
        "verify_simulink_release_encoding",
        ROOT / "tools" / "verify_simulink_release.py",
    )
    mex_dir = tmp_path / "mex"
    mex_dir.mkdir()
    for name in (
        *package_module.REQUIRED_MEX,
        *package_module.FULL_RUNTIME_DLLS,
    ):
        (mex_dir / name).write_bytes(fake_x64_pe())
    archive, _ = package_module.build_package(mex_dir, tmp_path / "out")
    matlab = tmp_path / "matlab.exe"
    matlab.write_bytes(b"MZ")

    class Result:
        returncode = 0
        stdout = "\u691c\u8a3c\u5b8c\u4e86 RADIA_IH_RELEASE_OK\n".encode("utf-8")
        stderr = b""

    monkeypatch.setattr(verify_module.subprocess, "run", lambda *args, **kwargs: Result())
    output = verify_module.run_matlab_smoke(archive, matlab)
    assert "RADIA_IH_RELEASE_OK" in output
    assert verify_module._console_safe("bad:\ufffd", "cp932") == "bad:\\ufffd"


def load_module(name, path):
    import importlib.util

    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


def fake_x64_pe():
    payload = bytearray(2048)
    payload[:2] = b"MZ"
    payload[60:64] = (128).to_bytes(4, "little")
    payload[128:132] = b"PE\0\0"
    payload[132:134] = (0x8664).to_bytes(2, "little")
    return bytes(payload)
