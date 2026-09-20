"""Exporter-owned GUI contracts; no Cubit or Qt process is launched."""
import importlib.util
import sys
import types
import pytest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / 'packages/cubit-mesh-export/src'
sys.path.insert(0, str(SOURCE))


def _load_from_source(module_path: str):
    """Load a module out of THIS checkout, whatever the install points at.

    conftest imports cubit_mesh_export before this module runs, so the package
    is already cached in sys.modules from wherever the editable install
    resolves -- a release worktree while the LAB pointer is drifted -- and the
    sys.path.insert above can no longer win.  The contract under test is about
    the repository layout, so the module is loaded by location instead of by
    name.  Editable drift is the drift checker's job, not this test's.
    """
    file = SOURCE / module_path
    spec = importlib.util.spec_from_file_location(
        "cubit_mesh_export_checkout_" + file.stem, file)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_bundled_fixture_is_product_owned_and_needs_no_radia():
    _find_sample_jou = _load_from_source(
        'cubit_mesh_export/smoke_test.py')._find_sample_jou
    path = _find_sample_jou()
    assert path.is_relative_to(SOURCE)
    text = path.read_text(encoding='utf-8')
    assert 'cubit-mesh-export' in text
    assert 'radia' not in text.lower()


def test_gui_startup_imports_only_local_menu(monkeypatch):
    path = SOURCE / 'cubit_mesh_export/cubit_gui/register_toolbar.py'
    calls = []
    monkeypatch.setitem(sys.modules, 'cubit_export_menu', types.SimpleNamespace(
        install_menu=lambda: calls.append(True) or True))
    original = sys.path[:]
    try:
        exec(compile(path.read_text(), str(path), 'exec'), {'__file__': str(path)})
    finally:
        sys.path[:] = original
    assert calls == [True]


def test_claro_actions_are_released_before_python_shutdown():
    source = (SOURCE / 'cubit_mesh_export/cubit_gui/cubit_export_menu.py').read_text(encoding='utf-8')
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
    startup = tmp_path / 'local/cubit-mesh-export/Cubit/startup.py'
    assert repr(sys.executable) in startup.read_text()
    assert 'cubit_gui/register_toolbar.py' in startup.read_text()
    startup.write_text('# stale startup')
    assert not installer.verify_panel_installation(verbose=False)[0]


def test_radia_contains_no_legacy_gui_bridges_or_test_implementation():
    for relative in ('install_panels.py', '_cubit_gui_compat.py',
                     'cubit_toolbar_smoke.py', 'panels/cubit_toolbar_probe.py',
                     'panels/startup.py', 'panels/register_toolbar.py',
                     'panels/cubit_export_menu.py'):
        assert not (ROOT / 'src/radia' / relative).exists(), relative
    assert not (ROOT / 'src/radia/panels/cubit_toolbar/toolbars/cubit_mesh_export_toolbar.ttb.tmpl').exists()


@pytest.mark.parametrize('name,fmt', [('netgen', 'netgen_vol'), ('gmsh', 'gmsh'),
                                    ('nastran', 'nastran'), ('vtk', 'vtk'),
                                    ('femeem', 'femeem'), ('meg', 'meg')])
@pytest.mark.parametrize('has_file', [False, True])
def test_official_toolbar_play_without_dunder_file(monkeypatch, tmp_path, name, fmt, has_file):
    source = SOURCE / f'cubit_mesh_export/cubit_gui/cubit_toolbar/scripts/export_{name}.py'
    (tmp_path / 'cubit_export_menu.py').write_text('# imported toolbar payload')
    monkeypatch.chdir(tmp_path)
    calls = []
    monkeypatch.setitem(sys.modules, 'cubit_export_menu', types.SimpleNamespace(launch_export=calls.append))
    namespace = {'__file__': str(tmp_path / source.name)} if has_file else {}
    original = sys.path[:]
    try:
        exec(compile(source.read_text(), str(source), 'exec'), namespace)
    finally:
        sys.path[:] = original
    assert calls == [fmt]
