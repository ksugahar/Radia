import hashlib
import importlib.util
import json
import sys
import tarfile
import xml.etree.ElementTree as ET
from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parents[1]
CME_SRC = PROJECT_ROOT / "packages" / "cubit-mesh-export" / "src"
if str(CME_SRC) not in sys.path:
    sys.path.insert(0, str(CME_SRC))


def _load_install_panels():
    path = PROJECT_ROOT / "src" / "radia" / "install_panels.py"
    spec = importlib.util.spec_from_file_location("radia_install_panels_test", path)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


def _fake_cubit(program_files: Path, version: str) -> Path:
    root = program_files / f"Coreform Cubit {version}"
    bin_dir = root / "bin"
    (bin_dir / "plugins").mkdir(parents=True)
    (bin_dir / "cubit.py").write_text("# fake cubit\n", encoding="utf-8")
    return root


def _patch_windows_env(monkeypatch, tmp_path: Path):
    monkeypatch.setenv("ProgramFiles", str(tmp_path / "ProgramFiles"))
    monkeypatch.setenv("ProgramFiles(x86)", "")
    monkeypatch.setenv("ProgramData", str(tmp_path / "ProgramData"))
    monkeypatch.setenv("LOCALAPPDATA", str(tmp_path / "LocalAppData"))
    monkeypatch.setenv("APPDATA", str(tmp_path / "AppData" / "Roaming"))
    monkeypatch.setenv("HOME", str(tmp_path / "Home"))
    monkeypatch.setenv("USERPROFILE", str(tmp_path / "Home"))
    monkeypatch.setenv("SystemDrive", str(tmp_path / "DriveC"))
    monkeypatch.delenv("CUBIT_PATH", raising=False)
    (tmp_path / "Home").mkdir()
    return Path(tmp_path / "ProgramFiles")


def test_find_cubit_bin_prefers_2025_12_over_2025_6(monkeypatch, tmp_path):
    install_panels = _load_install_panels()

    monkeypatch.setattr(install_panels.sys, "platform", "win32")
    program_files = _patch_windows_env(monkeypatch, tmp_path)
    _fake_cubit(program_files, "2025.3")
    _fake_cubit(program_files, "2025.6")
    expected = _fake_cubit(program_files, "2025.12") / "bin"

    assert Path(install_panels.find_cubit_bin()) == expected


def test_find_cubit_bin_rejects_pre_2025_12(monkeypatch, tmp_path):
    install_panels = _load_install_panels()

    monkeypatch.setattr(install_panels.sys, "platform", "win32")
    program_files = _patch_windows_env(monkeypatch, tmp_path)
    _fake_cubit(program_files, "2025.6")

    assert install_panels.find_cubit_bin() is None


def test_cubit_mesh_export_installer_requires_2025_12(monkeypatch, tmp_path):
    from cubit_mesh_export import install as cme_install

    monkeypatch.setattr(cme_install.sys, "platform", "win32")
    program_files = _patch_windows_env(monkeypatch, tmp_path)
    _fake_cubit(program_files, "2025.6")
    supported = _fake_cubit(program_files, "2025.12")

    assert cme_install._find_cubit_dir() == supported


def test_cubit_mesh_export_verify_requires_curver_pyd(monkeypatch, tmp_path):
    from cubit_mesh_export import install as cme_install

    pkg_dir = tmp_path / "pkg"
    cubit_dir = tmp_path / "Coreform Cubit 2025.12"
    plugins = cubit_dir / "bin" / "plugins"
    plugins.mkdir(parents=True)
    pkg_dir.mkdir()
    (pkg_dir / "cubit_mesh_export.ccm").write_bytes(b"ccm")
    nglib = tmp_path / "nglib.dll"
    ngcore = tmp_path / "ngcore.dll"
    nglib.write_bytes(b"nglib")
    ngcore.write_bytes(b"ngcore")
    monkeypatch.setattr(cme_install, "_find_netgen_dlls", lambda: (nglib, ngcore))

    ok, issues = cme_install.verify_deployment(pkg_dir, cubit_dir, verbose=False)

    assert not ok
    assert any("cubit_mesh_curver.pyd" in issue for issue in issues)


def test_helpers_only_install_does_not_touch_native_plugins(monkeypatch, tmp_path):
    from cubit_mesh_export import install as cme_install

    pkg_dir = tmp_path / "pkg"
    helpers_src = pkg_dir / "cubit_helpers"
    helpers_src.mkdir(parents=True)
    (helpers_src / "add_kelvin.py").write_text("VALUE = 1\n", encoding="utf-8")

    cubit_dir = tmp_path / "Coreform Cubit 2025.12"
    plugins = cubit_dir / "bin" / "plugins"
    plugins.mkdir(parents=True)
    native = plugins / "cubit_mesh_export.ccm"
    native.write_bytes(b"do-not-touch")

    monkeypatch.setattr(cme_install, "_package_dir", lambda: pkg_dir)
    monkeypatch.setattr(cme_install, "_find_cubit_dir", lambda: cubit_dir)
    monkeypatch.setattr(cme_install, "preflight", lambda *_args, **_kwargs: (True, []))

    assert cme_install.install_plugin(helpers_only=True) is True
    assert native.read_bytes() == b"do-not-touch"
    assert (
        plugins / "cubit_helpers" / "add_kelvin.py"
    ).read_text(encoding="utf-8") == "VALUE = 1\n"


def test_panel_startup_shim_is_generated_outside_package(monkeypatch, tmp_path):
    install_panels = _load_install_panels()

    monkeypatch.setattr(install_panels.sys, "platform", "win32")
    _patch_windows_env(monkeypatch, tmp_path)
    panels_dir = Path(install_panels._get_panels_dir())

    startup = Path(
        install_panels._generate_startup_script(str(panels_dir), all_users=True)
    )

    assert startup.parent == tmp_path / "ProgramData" / "Radia" / "Cubit"
    assert panels_dir not in startup.parents
    text = startup.read_text(encoding="utf-8")
    assert "register_toolbar.py" in text
    assert str(panels_dir).replace("\\", "/") in text


def test_install_panels_writes_and_verifies_current_user(monkeypatch, tmp_path):
    install_panels = _load_install_panels()

    monkeypatch.setattr(install_panels.sys, "platform", "win32")
    program_files = _patch_windows_env(monkeypatch, tmp_path)
    cubit_root = _fake_cubit(program_files, "2025.12")

    assert install_panels.install_panels(all_users=False) is True

    cubit_file = tmp_path / "Home" / ".cubit"
    startup = tmp_path / "LocalAppData" / "Radia" / "Cubit" / "radia_startup.py"
    toolbar_package = (
        tmp_path / "LocalAppData" / "Radia" / "Cubit"
        / "radia_export_toolbar.tar.gz"
    )
    ini = tmp_path / "AppData" / "Roaming" / "Coreform" / "Cubit.ini"

    assert cubit_file.is_file()
    assert startup.is_file()
    assert toolbar_package.is_file()
    assert "play" in cubit_file.read_text(encoding="utf-8")
    assert str(startup).replace("\\", "/") in cubit_file.read_text(encoding="utf-8")
    assert str(cubit_root / "bin" / "plugins").replace("\\", "/") in ini.read_text(
        encoding="utf-8"
    )
    ok, issues = install_panels.verify_panel_installation(all_users=False)
    assert ok, issues


def test_install_panels_refreshes_an_existing_official_toolbar(
        monkeypatch, tmp_path):
    install_panels = _load_install_panels()

    monkeypatch.setattr(install_panels.sys, "platform", "win32")
    program_files = _patch_windows_env(monkeypatch, tmp_path)
    _fake_cubit(program_files, "2025.12")
    toolbar_dir = (
        tmp_path / "LocalAppData" / "Radia" / "Cubit" / "Toolbars"
        / "radia_export_toolbar"
    )
    (toolbar_dir / "scripts").mkdir(parents=True)
    (toolbar_dir / "scripts" / "radia_export_menu.py").write_text(
        "# stale\n", encoding="utf-8")

    before = install_panels._verify_existing_toolbar_installations(
        all_users=False)
    assert any("stale" in issue for issue in before)
    assert install_panels.install_panels(all_users=False) is True

    issues = install_panels._verify_existing_toolbar_installations(
        all_users=False)
    assert issues == []
    menu = toolbar_dir / "scripts" / "radia_export_menu.py"
    assert menu.read_bytes() == (
        PROJECT_ROOT / "src" / "radia" / "panels"
        / "radia_export_menu.py"
    ).read_bytes()
    toolbar = toolbar_dir / "toolbars" / "radia_export_toolbar.ttb"
    assert toolbar_dir.as_posix() in toolbar.read_text(encoding="utf-8")


def test_official_toolbar_package_is_self_contained(tmp_path):
    install_panels = _load_install_panels()

    package = Path(install_panels.build_official_toolbar_package(tmp_path))
    assert package.name == "radia_export_toolbar.tar.gz"

    with tarfile.open(package, "r:gz") as archive:
        names = set(archive.getnames())
        assert "toolbars/radia_export_toolbar.ttb" in names
        assert "scripts/radia_export_menu.py" in names
        assert "icons/radia_export.svg" in names
        assert ".mappings" in names
        assert not any("__pycache__" in name or name.endswith(".pyc")
                       for name in names)
        toolbar_text = archive.extractfile(
            "toolbars/radia_export_toolbar.ttb").read().decode("utf-8")
        mappings = archive.extractfile(".mappings").read().decode("utf-8")

    root = ET.fromstring(toolbar_text)
    assert root.tag == "WorkflowToolbar"
    assert root.attrib == {
        "version": "1.0", "name": "Radia Export", "visible": "true",
    }
    buttons = root.findall("WTButton")
    assert len(buttons) == 6
    action_names = [
        button.find(".//WAction").attrib["name"] for button in buttons
    ]
    assert action_names == [
        "Netgen Vol (.vol)", "GMSH (.msh)", "Nastran (.bdf)",
        "VTK (.vtk)", "FEMEEM", "MEG (ELF/MAGIC)",
    ]

    referenced = [
        element.text for element in root.iter()
        if element.tag in {"filename", "icon"} and element.text
    ]
    assert all(f"{path} => " in mappings for path in referenced)


def test_native_build_is_worktree_relative_and_propagates_both_payloads():
    build_script = (
        PROJECT_ROOT / "src" / "cubit_plugin" / "cubit_build.ps1"
    ).read_text(encoding="utf-8")
    setup_script = (
        PROJECT_ROOT / "packages" / "cubit-mesh-export" / "setup.py"
    ).read_text(encoding="utf-8")

    assert r"S:\Radia\01_GitHub" not in build_script
    assert "$src = $PSScriptRoot" in build_script
    assert "packages\\cubit-mesh-export\\src\\cubit_mesh_export" in build_script
    for payload in ("cubit_mesh_export.ccm", "cubit_mesh_curver.pyd"):
        assert payload in build_script
        assert f'pkg_dir / "{payload}"' in setup_script


def test_distribution_ci_packages_the_exact_candidate_binaries():
    workflow = (
        PROJECT_ROOT / ".github" / "workflows" / "cubit-mesh-export.yml"
    ).read_text(encoding="utf-8")

    assert "CUBIT_MESH_EXPORT_SKIP_FRESHNESS_CHECK" not in workflow
    assert "python -m wheel tags" not in workflow
    assert "Root-Is-Purelib: false" in workflow
    assert "native_payloads.json" in workflow
    assert "download_release_asset.py" in workflow
    assert "Get-FileHash" in workflow
    assert '"netgen-mesher==6.2.2606"' in workflow
    assert '"ngsolve==6.2.2606"' in workflow
    assert "cubit_mesh_export/cubit_mesh_export.ccm" in workflow
    assert "cubit_mesh_export/cubit_mesh_curver.pyd" in workflow
    assert "tests\\test_cubit_menu_startup.py" in workflow
    assert "python tools\\audit_pyside6_only.py" in workflow

    setup_script = (
        PROJECT_ROOT / "packages" / "cubit-mesh-export" / "setup.py"
    ).read_text(encoding="utf-8")
    assert "class BinaryDistribution(Distribution)" in setup_script
    assert "setup(distclass=BinaryDistribution)" in setup_script

    package_dir = (
        PROJECT_ROOT / "packages" / "cubit-mesh-export" / "src"
        / "cubit_mesh_export"
    )
    manifest = json.loads(
        (package_dir / "native_payloads.json").read_text(encoding="utf-8"))
    payload = manifest["payloads"]["cubit_mesh_curver.pyd"]
    curver = (package_dir / "cubit_mesh_curver.pyd").read_bytes()
    assert payload["asset_name"].endswith(payload["sha256"] + ".pyd")
    assert payload["size"] == len(curver)
    assert payload["sha256"] == hashlib.sha256(curver).hexdigest()
