import hashlib
import ast
import importlib.util
import json
import os
import runpy
import shutil
import sys
import tarfile
import tomllib
import xml.etree.ElementTree as ET
from types import SimpleNamespace
from pathlib import Path

import pytest


PROJECT_ROOT = Path(__file__).resolve().parents[1]
CME_SRC = PROJECT_ROOT / "packages" / "cubit-mesh-export" / "src"
if str(CME_SRC) not in sys.path:
    sys.path.insert(0, str(CME_SRC))


@pytest.mark.parametrize("version,accepted", [("0.14.17", True), ("1.0.0", True),
                                               ("1.0.1", True), ("1.0.2", False)])
def test_radia_accepts_validated_exporter_versions_only(monkeypatch, version, accepted):
    import cubit_mesh_export
    from cubit_mesh_export import install
    tree = ast.parse((PROJECT_ROOT / "src/radia/__init__.py").read_text(encoding="utf-8"))
    names = {"__version__", "COMPAT_CUBIT_MESH_EXPORT_MIN", "COMPAT_CUBIT_MESH_EXPORT_MAX"}
    constants = {node.targets[0].id: ast.literal_eval(node.value)
                 for node in tree.body if isinstance(node, ast.Assign)
                 and isinstance(node.targets[0], ast.Name) and node.targets[0].id in names}
    monkeypatch.setitem(sys.modules, "radia", SimpleNamespace(**constants))
    monkeypatch.setattr(cubit_mesh_export, "__version__", version)
    assert install._check_radia_compat()[0] is accepted


@pytest.mark.parametrize("version,accepted", [("4.95.92", True), ("5.0.0", True),
                                              ("5.0.1", False), ("6.0.0", False)])
def test_exporter_radia5_compatibility_is_bounded(monkeypatch, version, accepted):
    import cubit_mesh_export
    from cubit_mesh_export import install

    monkeypatch.setitem(sys.modules, "radia", SimpleNamespace(
        __version__=version, COMPAT_CUBIT_MESH_EXPORT_MIN="0.5.0",
        COMPAT_CUBIT_MESH_EXPORT_MAX=cubit_mesh_export.__version__))
    assert install._check_radia_compat()[0] is accepted


def _load_install_panels():
    path = PROJECT_ROOT / "packages/cubit-mesh-export/src/cubit_mesh_export/toolbar_install.py"
    spec = importlib.util.spec_from_file_location("radia_install_panels_test", path)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


def _load_native_provenance():
    path = (
        PROJECT_ROOT / "packages" / "cubit-mesh-export" / "src"
        / "cubit_mesh_export" / "_native_provenance.py"
    )
    spec = importlib.util.spec_from_file_location(
        "cubit_mesh_export_native_provenance_test", path)
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

    assert startup.parent == tmp_path / "ProgramData" / "cubit-mesh-export" / "Cubit"
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
    startup = tmp_path / "LocalAppData" / "cubit-mesh-export" / "Cubit" / "startup.py"
    toolbar_package = (
        tmp_path / "LocalAppData" / "cubit-mesh-export" / "Cubit"
        / "cubit_mesh_export_toolbar.tar.gz"
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


def test_install_panels_deletes_only_exact_obsolete_exporter_assets(
        monkeypatch, tmp_path):
    install_panels = _load_install_panels()

    monkeypatch.setattr(install_panels.sys, "platform", "win32")
    program_files = _patch_windows_env(monkeypatch, tmp_path)
    _fake_cubit(program_files, "2025.12")

    old_startup_root = tmp_path / "LocalAppData" / "Radia" / "Cubit"
    old_toolbar = old_startup_root / "Toolbars" / "radia_export_toolbar"
    old_toolbar.mkdir(parents=True)
    (old_startup_root / "radia_startup.py").write_text(
        "obsolete\n", encoding="utf-8")
    (old_startup_root / "radia_export_toolbar.tar.gz").write_bytes(b"obsolete")
    (old_toolbar / "payload.py").write_text("obsolete\n", encoding="utf-8")
    sibling = old_startup_root / "keep-user-data.txt"
    sibling.write_text("keep\n", encoding="utf-8")
    old_settings = tmp_path / "AppData" / "Roaming" / "Radia"
    old_settings.mkdir(parents=True)
    (old_settings / "export_settings.json").write_text(
        "{}\n", encoding="utf-8")
    settings_sibling = old_settings / "keep-unrelated.json"
    settings_sibling.write_text("{}\n", encoding="utf-8")

    assert install_panels.install_panels(all_users=False) is True

    assert not (old_startup_root / "radia_startup.py").exists()
    assert not (old_startup_root / "radia_export_toolbar.tar.gz").exists()
    assert not old_toolbar.exists()
    assert not (old_settings / "export_settings.json").exists()
    assert sibling.read_text(encoding="utf-8") == "keep\n"
    assert settings_sibling.read_text(encoding="utf-8") == "{}\n"
    assert old_startup_root.is_dir()
    assert old_settings.is_dir()


def test_install_panels_refreshes_an_existing_official_toolbar(
        monkeypatch, tmp_path):
    install_panels = _load_install_panels()

    monkeypatch.setattr(install_panels.sys, "platform", "win32")
    program_files = _patch_windows_env(monkeypatch, tmp_path)
    _fake_cubit(program_files, "2025.12")
    toolbar_dir = (
        tmp_path / "LocalAppData" / "cubit-mesh-export" / "Cubit" / "Toolbars"
        / "cubit_mesh_export_toolbar"
    )
    (toolbar_dir / "scripts").mkdir(parents=True)
    (toolbar_dir / "scripts" / "cubit_export_menu.py").write_text(
        "# stale\n", encoding="utf-8")

    before = install_panels._verify_existing_toolbar_installations(
        all_users=False)
    assert any("stale" in issue for issue in before)
    assert install_panels.install_panels(all_users=False) is True

    issues = install_panels._verify_existing_toolbar_installations(
        all_users=False)
    assert issues == []
    menu = toolbar_dir / "scripts" / "cubit_export_menu.py"
    assert menu.read_bytes() == (
        PROJECT_ROOT / "packages/cubit-mesh-export/src/cubit_mesh_export/cubit_gui"
        / "cubit_export_menu.py"
    ).read_bytes()
    toolbar = toolbar_dir / "toolbars" / "cubit_mesh_export_toolbar.ttb"
    assert toolbar_dir.as_posix() in toolbar.read_text(encoding="utf-8")


def test_official_toolbar_package_is_self_contained(tmp_path):
    install_panels = _load_install_panels()

    package = Path(install_panels.build_official_toolbar_package(tmp_path))
    assert package.name == "cubit_mesh_export_toolbar.tar.gz"

    with tarfile.open(package, "r:gz") as archive:
        names = set(archive.getnames())
        assert "toolbars/cubit_mesh_export_toolbar.ttb" in names
        assert "scripts/cubit_export_menu.py" in names
        assert "icons/cubit_mesh_export.svg" in names
        assert ".mappings" in names
        assert not any("__pycache__" in name or name.endswith(".pyc")
                       for name in names)
        toolbar_text = archive.extractfile(
            "toolbars/cubit_mesh_export_toolbar.ttb").read().decode("utf-8")
        mappings = archive.extractfile(".mappings").read().decode("utf-8")

    root = ET.fromstring(toolbar_text)
    assert root.tag == "WorkflowToolbar"
    assert root.attrib == {
        "version": "1.0", "name": "Cubit Mesh Export", "visible": "true",
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
    assert all(path.startswith("./") for path in referenced)
    assert str(tmp_path).replace("\\", "/") not in toolbar_text
    assert str(tmp_path).replace("\\", "/") not in mappings


def test_unterminated_startup_block_never_discards_following_user_lines():
    install_panels = _load_install_panels()
    lines = [
        "set echo on\n",
        "## BEGIN cubit-mesh-export toolbar\n",
        'play "missing.py"\n',
        "user command that must survive\n",
    ]

    with pytest.raises(ValueError, match="unterminated"):
        install_panels._remove_existing_block(lines)

    assert lines[-1] == "user command that must survive\n"


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
    assert "_native_provenance.py" in build_script
    assert "record" in build_script
    assert "verify_manifest(repo_root, pkg_dir)" in setup_script
    assert "st_mtime" not in setup_script
    assert 'cmdclass={"build_py": CleanPackageBuild}' in setup_script


def test_package_build_deletes_stale_generated_modules(monkeypatch, tmp_path):
    import setuptools
    from setuptools.command.build_py import build_py

    package_root = tmp_path / "package"
    package_dir = package_root / "src" / "cubit_mesh_export"
    package_dir.mkdir(parents=True)
    original = PROJECT_ROOT / "packages" / "cubit-mesh-export"
    shutil.copy2(original / "setup.py", package_root / "setup.py")
    shutil.copy2(CME_SRC / "cubit_mesh_export" / "_native_provenance.py",
                 package_dir / "_native_provenance.py")
    monkeypatch.setenv("CUBIT_MESH_EXPORT_SKIP_FRESHNESS_CHECK", "1")
    setup_calls = []
    monkeypatch.setattr(setuptools, "setup", lambda **kwargs: setup_calls.append(kwargs))
    runpy.run_path(str(package_root / "setup.py"))

    generated = tmp_path / "build-lib" / "cubit_mesh_export"
    generated.mkdir(parents=True)
    stale = generated / "retired_module.py"
    stale.write_text("must disappear\n", encoding="utf-8")
    ran = []
    monkeypatch.setattr(build_py, "run", lambda _self: ran.append(True))
    command = setup_calls[0]["cmdclass"]["build_py"](setuptools.Distribution())
    command.build_lib = str(tmp_path / "build-lib")
    command.run()

    assert ran == [True]
    assert not generated.exists()


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
    # PowerShell continuation and argument order do not define the contract.
    install_commands = [line.split("pip install", 1)[1].split()
                        for line in workflow.replace("`\n", " ").splitlines()
                        if "python -m pip install " in line]
    required_test_dependencies = {"setuptools", "pytest", "numpy", "psutil"}
    assert any(required_test_dependencies <= set(args) for args in install_commands)
    assert "_native_provenance.py" in workflow
    assert "verify --repo-root . --package-dir $destination" in workflow
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
    assert "setup(distclass=BinaryDistribution," in setup_script
    assert 'cmdclass={"build_py": CleanPackageBuild}' in setup_script

    package_dir = (
        PROJECT_ROOT / "packages" / "cubit-mesh-export" / "src"
        / "cubit_mesh_export"
    )
    manifest = json.loads(
        (package_dir / "native_payloads.json").read_text(encoding="utf-8"))
    provenance = _load_native_provenance()
    assert manifest["schema"] == provenance.SCHEMA
    source_hash, source_count = provenance.native_source_digest(PROJECT_ROOT)
    assert manifest["source"]["hash_format"] == provenance.HASH_FORMAT
    assert manifest["source"]["tree_sha256"] == source_hash
    assert manifest["source"]["file_count"] == source_count
    assert provenance.verify_manifest(PROJECT_ROOT, package_dir) == []
    for name in provenance.REQUIRED_PAYLOADS:
        payload = manifest["payloads"][name]
        content = (package_dir / name).read_bytes()
        assert payload["size"] == len(content)
        assert payload["sha256"] == hashlib.sha256(content).hexdigest()
    curver = manifest["payloads"]["cubit_mesh_curver.pyd"]
    assert curver["asset_name"].endswith(curver["sha256"] + ".pyd")


def test_native_source_provenance_is_content_addressed_not_timestamped(tmp_path):
    provenance = _load_native_provenance()
    source = tmp_path / "src" / "cubit_plugin" / "native.cpp"
    source.parent.mkdir(parents=True)
    source.write_text("int answer = 42;\n", encoding="utf-8")

    before = provenance.native_source_digest(tmp_path)
    stat = source.stat()
    os.utime(source, (stat.st_atime + 100, stat.st_mtime + 100))
    assert provenance.native_source_digest(tmp_path) == before

    source.write_bytes(b"int answer = 42;\r\n")
    assert provenance.native_source_digest(tmp_path) == before

    source.write_text("int answer = 43;\n", encoding="utf-8")
    assert provenance.native_source_digest(tmp_path) != before


def test_native_manifest_rejects_drift_in_either_required_payload(
        monkeypatch, tmp_path):
    provenance = _load_native_provenance()
    source = tmp_path / "src" / "cubit_plugin" / "native.cpp"
    source.parent.mkdir(parents=True)
    source.write_text("int answer = 42;\n", encoding="utf-8")
    package_dir = tmp_path / "package"
    package_dir.mkdir()
    for name in provenance.REQUIRED_PAYLOADS:
        (package_dir / name).write_bytes((name + "\n").encode())
    (package_dir / "native_payloads.json").write_text(
        '{"payloads": {}}', encoding="utf-8")
    monkeypatch.setattr(provenance, "_source_commit", lambda _root: "a" * 40)
    provenance.record_manifest(tmp_path, package_dir)
    assert provenance.verify_manifest(tmp_path, package_dir) == []

    for name in provenance.REQUIRED_PAYLOADS:
        path = package_dir / name
        original = path.read_bytes()
        path.write_bytes(original + b"drift")
        errors = provenance.verify_manifest(tmp_path, package_dir)
        assert any(name in error for error in errors)
        path.write_bytes(original)


@pytest.mark.parametrize("damage", [
    None, "schema", "ccm-content", "pyd-content", "ccm-missing", "pyd-missing",
    "size",
])
def test_sdist_setup_verifies_payloads_without_native_sources(
    monkeypatch, tmp_path, damage
):
    import setuptools

    provenance = _load_native_provenance()
    package_root = tmp_path / "unpacked" / "cubit-mesh-export"
    package_dir = package_root / "src" / "cubit_mesh_export"
    package_dir.mkdir(parents=True)
    for name in provenance.REQUIRED_PAYLOADS:
        (package_dir / name).write_bytes(name.encode())
    manifest_path = package_dir / "native_payloads.json"
    manifest_path.write_text('{"payloads": {}}', encoding="utf-8")
    monkeypatch.setattr(provenance, "_source_commit", lambda _root: "a" * 40)
    provenance.record_manifest(tmp_path, package_dir)
    assert provenance.verify_manifest(None, package_dir) == []

    if damage in ("schema", "size"):
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
        if damage == "schema":
            manifest["schema"] = "unsupported"
        else:
            manifest["payloads"]["cubit_mesh_export.ccm"]["size"] += 1
        manifest_path.write_text(json.dumps(manifest), encoding="utf-8")
    elif damage is not None:
        name = ("cubit_mesh_export.ccm" if damage.startswith("ccm")
                else "cubit_mesh_curver.pyd")
        path = package_dir / name
        if damage.endswith("missing"):
            path.unlink()
        else:
            # Same size: the digest must catch corruption independently.
            data = path.read_bytes()
            path.write_bytes(bytes([data[0] ^ 1]) + data[1:])

    original_package = PROJECT_ROOT / "packages" / "cubit-mesh-export"
    shutil.copy2(original_package / "setup.py", package_root / "setup.py")
    shutil.copy2(CME_SRC / "cubit_mesh_export" / "_native_provenance.py",
                 package_dir / "_native_provenance.py")
    assert not (tmp_path / "src" / "cubit_plugin").exists()
    monkeypatch.delenv("CUBIT_MESH_EXPORT_SKIP_FRESHNESS_CHECK", raising=False)
    setup_calls = []
    monkeypatch.setattr(setuptools, "setup", lambda **kwargs: setup_calls.append(kwargs))
    if damage is None:
        runpy.run_path(str(package_root / "setup.py"))
        assert len(setup_calls) == 1
    else:
        with pytest.raises(SystemExit) as error:
            runpy.run_path(str(package_root / "setup.py"))
        assert error.value.code == 1
        assert not setup_calls


def test_distribution_metadata_matches_the_only_supported_native_wheel():
    pyproject = tomllib.loads((
        PROJECT_ROOT / "packages" / "cubit-mesh-export" / "pyproject.toml"
    ).read_text(encoding="utf-8"))
    project = pyproject["project"]
    assert project["requires-python"] == ">=3.12,<3.13"
    assert "Programming Language :: Python :: 3.12" in project["classifiers"]
    assert "Operating System :: Microsoft :: Windows" in project["classifiers"]
