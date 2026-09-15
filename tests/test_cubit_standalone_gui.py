"""Exporter-owned GUI contracts; no Cubit or Qt process is launched."""
import importlib.util
import sys
import types
import runpy
import pytest
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


def test_legacy_bridge_loads_assets_without_importing_native_packages(monkeypatch, tmp_path):
    bridge = runpy.run_path(str(ROOT / 'src/radia/_cubit_gui_compat.py'))
    root = tmp_path / 'cubit_mesh_export'
    root.mkdir()
    (root / 'toolbar_install.py').write_text('VALUE = 42\n', encoding='utf-8')
    distribution = types.SimpleNamespace(locate_file=lambda name: tmp_path / name,
                                         read_text=lambda name: None)
    monkeypatch.setattr(bridge['importlib'].metadata, 'distribution', lambda name: distribution)
    namespace = {'__name__': 'legacy_test'}
    bridge['load_exporter_file']('toolbar_install.py', namespace)
    assert namespace['VALUE'] == 42
    assert Path(namespace['__file__']) == root / 'toolbar_install.py'


def test_legacy_cubit_without_external_python_fails_with_reinstall_guidance(monkeypatch):
    bridge = runpy.run_path(str(ROOT / 'src/radia/_cubit_gui_compat.py'))
    def absent(name):
        raise bridge['importlib'].metadata.PackageNotFoundError(name)
    monkeypatch.setattr(bridge['importlib'].metadata, 'distribution', absent)
    monkeypatch.delenv('CUBIT_MESH_EXPORT_PYTHON', raising=False)
    monkeypatch.delenv('RADIA_PYTHON', raising=False)
    with pytest.raises(ModuleNotFoundError, match='run cubit-plugin-install'):
        bridge['load_exporter_file']('cubit_gui/register_toolbar.py', {})


def test_radia_paths_are_bridges_not_duplicate_gui_implementations():
    for relative in ('install_panels.py', 'panels/register_toolbar.py', 'panels/radia_export_menu.py'):
        text = (ROOT / 'src/radia' / relative).read_text(encoding='utf-8')
        assert 'load_exporter_file' in text
        assert len(text.splitlines()) < 10
    assert not (ROOT / 'src/radia/panels/cubit_toolbar/toolbars/radia_export_toolbar.ttb.tmpl').exists()


@pytest.mark.parametrize('name,fmt', [('netgen', 'netgen_vol'), ('gmsh', 'gmsh'),
                                    ('nastran', 'nastran'), ('vtk', 'vtk'),
                                    ('femeem', 'femeem'), ('meg', 'meg')])
@pytest.mark.parametrize('has_file', [False, True])
def test_official_toolbar_play_without_dunder_file(monkeypatch, tmp_path, name, fmt, has_file):
    source = SOURCE / f'cubit_mesh_export/cubit_gui/cubit_toolbar/scripts/export_{name}.py'
    (tmp_path / 'radia_export_menu.py').write_text('# imported toolbar payload')
    monkeypatch.chdir(tmp_path)
    calls = []
    monkeypatch.setitem(sys.modules, 'radia_export_menu', types.SimpleNamespace(launch_export=calls.append))
    namespace = {'__file__': str(tmp_path / source.name)} if has_file else {}
    original = sys.path[:]
    try:
        exec(compile(source.read_text(), str(source), 'exec'), namespace)
    finally:
        sys.path[:] = original
    assert calls == [fmt]
