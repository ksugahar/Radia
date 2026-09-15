"""Exporter-owned GUI contracts; no Cubit or Qt process is launched."""
import importlib.util
import sys
import types
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / 'packages/cubit-mesh-export/src'
sys.path.insert(0, str(SOURCE))


def test_bundled_fixture_needs_no_radia():
    from cubit_mesh_export.smoke_test import _find_sample_jou
    path = _find_sample_jou()
    assert path.is_relative_to(SOURCE)
    assert path.read_text(encoding='utf-8') == (ROOT / 'src/radia/panels/samples/ih_bem_sample.jou').read_text(encoding='utf-8')


def test_gui_startup_imports_only_local_menu(monkeypatch):
    path = SOURCE / 'cubit_mesh_export/cubit_gui/register_toolbar.py'
    calls = []
    monkeypatch.setitem(sys.modules, 'radia_export_menu', types.SimpleNamespace(
        install_menu=lambda: calls.append(True) or True))
    original = sys.path[:]
    try:
        exec(compile(path.read_text(), str(path), 'exec'), {'__file__': str(path)})
    finally:
        sys.path[:] = original
    assert calls == [True]


def test_claro_actions_are_released_before_python_shutdown():
    source = (SOURCE / 'cubit_mesh_export/cubit_gui/radia_export_menu.py').read_text(encoding='utf-8')
    assert 'app.aboutToQuit.connect(remove_menu)' in source
    assert 'emclaro.remove_menu_items(_CLARO_COMPONENT)' in source


def test_standalone_registration_in_scratch_profile(monkeypatch, tmp_path):
    from cubit_mesh_export import toolbar_install as installer
    cubit = tmp_path / 'Coreform Cubit 2025.12/bin'
    cubit.mkdir(parents=True)
    (cubit / 'cubit.py').write_text('')
    monkeypatch.setenv('CUBIT_PATH', str(cubit))
    monkeypatch.setenv('LOCALAPPDATA', str(tmp_path / 'local'))
    monkeypatch.setenv('APPDATA', str(tmp_path / 'roaming'))
    monkeypatch.setenv('USERPROFILE', str(tmp_path / 'home'))
    monkeypatch.setenv('HOME', str(tmp_path / 'home'))
    assert installer.install_panels()
    assert installer.verify_panel_installation(verbose=False) == (True, [])
    startup = tmp_path / 'local/Radia/Cubit/radia_startup.py'
    assert repr(sys.executable) in startup.read_text()
    assert 'cubit_gui/register_toolbar.py' in startup.read_text()
    startup.write_text('# stale startup')
    assert not installer.verify_panel_installation(verbose=False)[0]
