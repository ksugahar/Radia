from pathlib import Path
import zipfile

import pytest


ROOT = Path(__file__).resolve().parents[1]


def test_package_builder_requires_native_ih_assets():
    text = (ROOT / "tools" / "package_simulink_release.py").read_text(encoding="utf-8")
    assert "radia_ih.slx" in text
    assert "radia_ih_eddy_sfun.mexw64" in text
    assert "radia_ih_thermal_sfun.mexw64" in text
    assert "validateIHNativeConfig.m" in text
    assert "manifest.json" in text
    assert "SHA256SUMS.txt" in text


def test_package_builder_fails_when_mex_is_missing(tmp_path):
    import importlib.util

    path = ROOT / "tools" / "package_simulink_release.py"
    spec = importlib.util.spec_from_file_location("package_simulink_release", path)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    with pytest.raises(FileNotFoundError, match="Required native IH binary"):
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
    for name in package_module.REQUIRED_MEX:
        (mex_dir / name).write_bytes(fake_x64_pe())
    archive, sums = package_module.build_package(mex_dir, tmp_path / "out")
    manifest = verify_module.verify_archive(archive)
    assert sums.read_text(encoding="ascii").split()[0] == \
        package_module.sha256(archive)
    assert manifest["release_channel"] == "preview"
    assert manifest["operator_assembly"] == "preassembled"
    assert (tmp_path / "out" / "manifest.json").is_file()
    assert len(sums.read_text(encoding="ascii").splitlines()) == 2
    with zipfile.ZipFile(archive) as bundle:
        names = set(bundle.namelist())
    assert "matlab/radia_ih.slx" in names
    assert not any("LUT" in name or "makeIHPlant" in name for name in names)
    assert "matlab/radia_simulink_library.slx" not in names


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
