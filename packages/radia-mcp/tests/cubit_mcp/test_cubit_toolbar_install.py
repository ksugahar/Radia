"""Guards: installer, toolbar installer and MCP session agree on the supported Cubit
release; the toolbar payload is self-contained (menu, Kelvin helper, launcher, icon);
the startup shim resolves radia_mcp's GUI through RADIA_MCP_CUBIT_PYTHON; and a
reinstall replaces its own startup block and Cubit.ini toolbar entry.

Assumes: a fake Cubit install tree and profile under tmp_path; no Cubit runs and
no real user profile is touched (HOME/APPDATA/LOCALAPPDATA point into tmp_path).
"""
from __future__ import annotations

import ast
import sys
import tarfile
import types
from pathlib import Path

import pytest

from radia_mcp.cubit import install, layout, session, toolbar_install

CUBIT = Path(layout.__file__).parent

# (directory name, BUILD_INFO text or None, supported)
RELEASE_CASES = [
    ("Coreform Cubit 2026.8", None, False),              # name only
    ("Coreform Cubit 2026.8", "short=2026.6\n", False),  # name lies
    ("renamed-cubit", "short=2026.8\r\n", True),        # renamed, exact release
    ("renamed-cubit", "# short=2026.8\nshort = 2026.8 \n", True),
    ("renamed-cubit", "version=2026.8\n", False),       # no short field
]


def _install_root(tmp_path, name, build_info):
    root = tmp_path / name
    (root / "bin" / "plugins").mkdir(parents=True)
    (root / "plugins").mkdir()
    (root / "bin" / "cubitx.exe").touch()
    (root / "bin" / "cubit.py").touch()
    if build_info is not None:
        (root / "BUILD_INFO").write_bytes(build_info.encode("utf-8"))
    return root


@pytest.mark.parametrize("init_present,native_present,expected", [
    (False, False, False), (True, False, False), (False, True, False), (True, True, True),
])
def test_packaged_binding_requires_both_files(tmp_path, init_present, native_present, expected):
    package = tmp_path / "cubit"
    package.mkdir()
    if init_present:
        (package / "__init__.py").touch()
    if native_present:
        (package / "_cubit3.pyd").touch()
    assert layout.has_python_binding(tmp_path) is expected


@pytest.mark.parametrize("name,build_info,supported", RELEASE_CASES)
def test_all_entry_points_agree_on_the_release(tmp_path, monkeypatch, name, build_info, supported):
    root = _install_root(tmp_path, name, build_info)
    monkeypatch.delenv("CUBIT_INSTALL_DIR", raising=False)
    monkeypatch.setenv("CUBIT_BIN_DIR", str(root / "bin"))
    # Ignore a real Cubit on this machine: no install globs, no registry keys.
    monkeypatch.setattr(session.glob, "glob", lambda pattern: [])
    registry = types.ModuleType("winreg")
    registry.HKEY_LOCAL_MACHINE = 0

    def no_key(*args):
        raise OSError("no registry in this test")

    registry.OpenKey = no_key
    monkeypatch.setitem(sys.modules, "winreg", registry)
    verdicts = {
        "layout": layout.is_supported_cubit_root(root),
        "installer": install._find_cubit_dir() == root.resolve(),
        "toolbar": toolbar_install.find_cubit_bin() == (str((root / "bin").resolve()) if supported else None),
        "session": session.find_cubit_install(str(root / "bin")) == (root / "bin").resolve(),
    }
    expected = {"layout": supported, "installer": supported, "toolbar": True, "session": supported}
    assert verdicts == expected, (name, build_info, verdicts)


def test_plugin_is_deployed_under_the_loader_directory(tmp_path):
    root = _install_root(tmp_path, "renamed-cubit", "short=2026.8\n")
    assert layout.command_plugin_path(root / "bin") == root / "plugins" / "cubit_to_vol_export.dll"
    assert layout.COMMAND_PLUGIN == install.COMMAND_PLUGIN


def test_ported_modules_never_read_build_info_themselves():
    for name in ("install.py", "toolbar_install.py", "toolbar_smoke.py", "smoke_test.py"):
        tree = ast.parse((CUBIT / name).read_text(encoding="utf-8"))
        assert not any(isinstance(n, ast.Constant) and n.value == "BUILD_INFO"
                       for n in ast.walk(tree)), name


def test_install_requires_a_built_plugin_directory(tmp_path, monkeypatch):
    monkeypatch.delenv("RADIA_CUBIT_PLUGIN_DIR", raising=False)
    with pytest.raises(RuntimeError, match="plugin directory not set"):
        install.resolve_plugin_dir(None)
    with pytest.raises(RuntimeError, match="cubit_to_vol_export.dll not found"):
        install.resolve_plugin_dir(tmp_path)
    (tmp_path / "cubit_to_vol_export.dll").write_bytes(b"MZ")
    assert install.resolve_plugin_dir(tmp_path) == tmp_path.resolve()


def test_toolbar_package_is_self_contained(tmp_path):
    package = Path(toolbar_install.build_official_toolbar_package(output_dir=tmp_path))
    with tarfile.open(package) as archive:
        names = sorted(m.name for m in archive.getmembers() if m.isfile())
        mappings = archive.extractfile(".mappings").read().decode()
        toolbar = archive.extractfile("toolbars/vol_export_toolbar.ttb").read().decode()
    assert names == [".mappings", "icons/vol_export.svg", "scripts/cubit_export_menu.py",
                     "scripts/export_vol.py", "scripts/kelvin.py",
                     "toolbars/vol_export_toolbar.ttb"], names
    assert "./scripts/kelvin.py => scripts/kelvin.py" in mappings, mappings
    assert "@TOOLBAR_INSTALL_DIR@" not in toolbar


@pytest.fixture
def profile(tmp_path, monkeypatch):
    home = tmp_path / "home"
    for name, value in (("LOCALAPPDATA", home / "AppData/Local"),
                        ("APPDATA", home / "AppData/Roaming")):
        value.mkdir(parents=True)
        monkeypatch.setenv(name, str(value))
    monkeypatch.setenv("USERPROFILE", str(home))
    monkeypatch.setenv("HOME", str(home))
    return home


def test_startup_shim_resolves_the_gui_through_the_recorded_python(profile):
    path = Path(toolbar_install._generate_startup_script(toolbar_install._get_panels_dir()))
    expected = (profile / "AppData/Local/radia-mcp/Cubit" if sys.platform == "win32"
                else profile / ".radia-mcp/cubit")
    assert path.parent == expected, (path, expected)
    text = path.read_text(encoding="utf-8")
    assert "RADIA_MCP_CUBIT_PYTHON" in text
    register = (CUBIT / "gui" / "register_toolbar.py").as_posix()
    assert register in text.replace("\\\\", "/"), text
    # The resolver probe, run by the recorded interpreter, finds this checkout's GUI.
    body = next(ast.literal_eval(n.value.args[0]) for n in ast.parse(text).body[-1:])
    probe = ast.parse(body)
    assert any(isinstance(n, ast.Constant) and "find_spec" in str(n.value) for n in ast.walk(probe))


def test_owned_startup_block_and_ini_entry_are_replaced(profile):
    cubit_file = profile / ".cubit"
    cubit_file.write_bytes(b"journal off\n## BEGIN radia-mcp cubit toolbar\nplay \"x\"\n"
                           b"## END radia-mcp cubit toolbar\nkeep me\n")
    lines, _, _ = toolbar_install._read_startup_for_edit(cubit_file)
    assert "".join(lines) == "journal off\nkeep me\n", lines
    ini = profile / "AppData/Roaming/Coreform/Cubit.ini"
    ini.parent.mkdir(parents=True)
    ini.write_text("[workflow]\nconverted_old=true\n"
                   "toolbars\\0=C:/other/toolbars/other.ttb\n"
                   "toolbars\\1=C:/u/AppData/Local/radia-mcp/Cubit/Toolbars/"
                   "vol_export_toolbar/toolbars/vol_export_toolbar.ttb\n",
                   encoding="utf-8")
    assert toolbar_install._register_official_toolbar(str(ini)) is True
    _, entries = toolbar_install._workflow_toolbar_entries(ini.read_text(encoding="utf-8"))
    values = [value for _, value, _ in entries]
    target = (toolbar_install._profile_toolbar_dir(str(ini)) / "toolbars"
              / "vol_export_toolbar.ttb").as_posix()
    assert values == ["C:/other/toolbars/other.ttb", target], values
    assert toolbar_install._ini_has_official_toolbar(str(ini))
    assert toolbar_install._register_official_toolbar(str(ini)) is False  # idempotent


@pytest.mark.parametrize("argv", [["--help"], ["--verify"], ["--uninstall", "--verify-only"]])
def test_unknown_or_conflicting_flags_never_install(monkeypatch, argv):
    calls = []
    for name in ("install_panels", "uninstall_panels", "verify_panel_installation",
                 "build_official_toolbar_package"):
        monkeypatch.setattr(toolbar_install, name, lambda *a, _n=name, **k: calls.append(_n))
    with pytest.raises(SystemExit):
        toolbar_install.main(argv)
    assert calls == [], calls


def test_export_plugin_directory_prefers_the_build_and_needs_the_plugin(tmp_path, monkeypatch):
    """Scripts pass this directory to -commandplugindir: a build.py output when
    RADIA_CUBIT_PLUGIN_DIR is set, else the installed <root>/plugins (not bin/plugins)."""
    root = tmp_path / "cubit"
    (root / "bin" / "plugins").mkdir(parents=True)
    (root / "plugins").mkdir()
    (root / "bin" / "cubitx.exe").touch()
    built = tmp_path / "built"
    built.mkdir()
    monkeypatch.delenv("RADIA_CUBIT_PLUGIN_DIR", raising=False)
    with pytest.raises(FileNotFoundError, match="cubit_to_vol_export.dll"):
        layout.export_plugin_directory(root / "bin")
    (root / "bin" / "plugins" / layout.COMMAND_PLUGIN).touch()   # the wrong place
    with pytest.raises(FileNotFoundError):
        layout.export_plugin_directory(root / "bin")
    (root / "plugins" / layout.COMMAND_PLUGIN).touch()
    assert layout.export_plugin_directory(root / "bin") == root / "plugins"
    monkeypatch.setenv("RADIA_CUBIT_PLUGIN_DIR", str(built))
    with pytest.raises(FileNotFoundError):
        layout.export_plugin_directory(root / "bin")
    (built / layout.COMMAND_PLUGIN).touch()
    assert layout.export_plugin_directory(root / "bin") == built.resolve()
