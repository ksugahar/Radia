"""Regression tests for Cubit's startup hook and GUI ownership boundary."""

from __future__ import annotations

import ast
import importlib.util
import os
from pathlib import Path
import sys
import types

import pytest

ROOT = Path(__file__).resolve().parents[1]
REGISTER_TOOLBAR = ROOT / "packages/cubit-mesh-export/src/cubit_mesh_export/cubit_gui/register_toolbar.py"
EXPORT_MENU = ROOT / "packages/cubit-mesh-export/src/cubit_mesh_export/cubit_gui/cubit_export_menu.py"
AUDIT_SCRIPT = ROOT / "tools" / "audit_pyside6_only.py"

_AUDIT_SPEC = importlib.util.spec_from_file_location(
    "cubit_export_audit_pyside6_only", AUDIT_SCRIPT
)
assert _AUDIT_SPEC is not None and _AUDIT_SPEC.loader is not None
_AUDIT_MODULE = importlib.util.module_from_spec(_AUDIT_SPEC)
_AUDIT_SPEC.loader.exec_module(_AUDIT_MODULE)
check_deployed_panel_source = _AUDIT_MODULE.check_deployed_panel_source
check_pyside6_ownership_boundary = _AUDIT_MODULE.check_pyside6_ownership_boundary


def test_pyside6_is_confined_to_cubit_owned_code_and_tests():
    assert check_pyside6_ownership_boundary() == []


@pytest.mark.parametrize("check", [
    _AUDIT_MODULE.check_no_legacy_qt_imports,
    check_pyside6_ownership_boundary,
])
@pytest.mark.parametrize("returncode", [0, 1, 128])
def test_import_audit_distinguishes_no_matches_from_git_failure(monkeypatch, check, returncode):
    monkeypatch.setattr(_AUDIT_MODULE.subprocess, "run", lambda *a, **kw:
                        types.SimpleNamespace(returncode=returncode, stdout="",
                                              stderr="repository unavailable"))
    issues = check()
    if returncode == 128:
        assert issues == ["git audit failed (128): repository unavailable"]
    else:
        assert issues == []


@pytest.mark.parametrize("returncode", [2, -1073741819])
def test_failed_export_process_cannot_pass_with_an_existing_artifact(monkeypatch, tmp_path, returncode):
    path = EXPORT_MENU.parents[1] / "smoke_test.py"
    spec = importlib.util.spec_from_file_location("cubit_process_smoke", path)
    smoke = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(smoke)
    monkeypatch.setattr(smoke, "_find_cubit_exe", lambda: tmp_path / "coreform_cubit.exe")
    monkeypatch.setattr(smoke, "_find_sample_jou", lambda *_: tmp_path / "sample.jou")
    monkeypatch.setattr(smoke.tempfile, "mkdtemp", lambda **_: str(tmp_path))

    def failed_export(*args, **kwargs):
        (tmp_path / "smoke.vol").write_text("exported before failure", encoding="utf-8")
        return types.SimpleNamespace(returncode=returncode, stdout="export finished", stderr="shutdown failed")

    monkeypatch.setattr(smoke.subprocess, "run", failed_export)
    monkeypatch.setattr(smoke, "_validate_exported_vol", lambda *a, **kw:
                        pytest.fail("failed process must not be accepted via artifact validation"))
    assert smoke.run_smoke_test() == 1
    assert (tmp_path / "smoke.vol").is_file()
    assert "shutdown failed" in (tmp_path / "cubit.log").read_text(encoding="utf-8")


def test_deployment_audit_rejects_startup_from_another_checkout(tmp_path):
    gui_rel = Path("packages/cubit-mesh-export/src/cubit_mesh_export/cubit_gui")
    expected = tmp_path / "current" / gui_rel / "register_toolbar.py"
    foreign = tmp_path / "old-release" / gui_rel / "register_toolbar.py"
    expected.parent.mkdir(parents=True)
    foreign.parent.mkdir(parents=True)
    expected.write_text("# current\n", encoding="utf-8")
    foreign.write_text("# stale\n", encoding="utf-8")

    startup = tmp_path / "startup.py"
    startup.write_text(
        f"exec(open(r'{foreign.as_posix()}').read())\n",
        encoding="utf-8",
    )
    cubit_file = tmp_path / ".cubit"
    cubit_file.write_text(
        "## BEGIN cubit-mesh-export toolbar\n"
        f'play "{startup.as_posix()}"\n'
        "## END cubit-mesh-export toolbar\n",
        encoding="utf-8",
    )

    status, issues = check_deployed_panel_source(cubit_file, expected)
    assert status == "checked"
    assert len(issues) == 1
    assert "different checkout" in issues[0]
    assert expected.as_posix().lower() in issues[0]

    startup.write_text(
        f"exec(open(r'{expected.as_posix()}').read())\n",
        encoding="utf-8",
    )
    status, issues = check_deployed_panel_source(cubit_file, expected)
    assert status == "checked"
    assert issues == []


def test_cubit_startup_installs_menu_through_claro_api():
    """The startup hook must install the menu via Cubit's Claro API.

    History: the menu used to be injected straight into Cubit's QMenuBar
    with PySide6.  Cubit discarded it during cold start because it did
    not own it, so the menu "disappeared".  MEASURED on Cubit 2025.12:
    ``emclaro.is_loaded()`` is already True inside ``~/.cubit``, and a
    menu registered through ``emclaro.add_to_menu()`` is still present
    once startup completes.
    """
    source = REGISTER_TOOLBAR.read_text(encoding="utf-8")

    assert "cubit_export_menu.install_menu()" in source
    assert "import radia" not in source.replace("import cubit_export_menu", "")


def test_export_menu_uses_claro_api_not_qmenubar_injection():
    """The menu module must register through emclaro, not QMenuBar."""
    source = EXPORT_MENU.read_text(encoding="utf-8")

    assert "import emclaro" in source
    assert "add_to_menu" in source
    assert "remove_menu_items" in source
    # No QMenuBar injection: that is the bug this replaced.
    assert "menu_bar.addMenu" not in source
    assert "QMenu(" not in source


def test_export_menu_does_not_force_a_journal_save():
    """A loaded Cubit model must export without creating a root-side log."""
    source = EXPORT_MENU.read_text(encoding="utf-8")

    assert "ensure_jou_path" not in source
    assert "getSaveFileName" not in source
    assert "save journal" not in source
    assert "_current_journal_hint" in source


def test_claro_export_menu_not_removed_through_qt():
    """"Cubit Mesh Export" must not be torn down via Qt.

    Enumerating Claro-owned QMenu/QAction objects from PySide6 deadlocks
    Cubit (reproduced twice: the GUI stops responding and the play script
    never returns).  Cleanup goes through emclaro.remove_menu_items().
    """
    source = REGISTER_TOOLBAR.read_text(encoding="utf-8")

    assert '"Export Mesh", "Cubit Mesh Export"' not in source
    assert '"Cubit Mesh Export",' not in source


def test_find_claro_matches_capital_c_object_name():
    """Cubit 2025.12 names its main window ``Claro`` -- capital C.

    The previous lowercase comparison never matched, and the
    any-QMainWindow fallback then returned ``QtJournalEditor`` instead,
    so dialogs were parented to the wrong window and legacy-menu cleanup
    operated on the wrong menu bar.
    """
    for path in (EXPORT_MENU, REGISTER_TOOLBAR):
        source = path.read_text(encoding="utf-8")
        if "objectName()" not in source:
            continue
        assert 'objectName() == "claro"' not in source, path
        assert 'objectName().lower() == "claro"' in source, path
        assert "QtJournalEditor" in source, path


def test_nastran_action_uses_the_solver_neutral_command():
    """The visible menu uses the package-owned BDF command."""
    source = EXPORT_MENU.read_text(encoding="utf-8")

    assert "export nastran_bdf" in source


def test_claro_activation_strings_dispatch_every_export_format(monkeypatch):
    """Cubit must be able to execute every PyAction activation string.

    ``setActivateMethod`` stores source text for Cubit's embedded Python,
    rather than a Python callable.  Extract only the source generator so this
    contract stays testable without importing Cubit's private PySide6 runtime.
    """
    source = EXPORT_MENU.read_text(encoding="utf-8")
    tree = ast.parse(source, filename=str(EXPORT_MENU))
    activate_node = next(
        node for node in tree.body
        if isinstance(node, ast.FunctionDef) and node.name == "_activate_code"
    )
    namespace = {"os": os, "__file__": str(EXPORT_MENU)}
    exec(compile(ast.Module([activate_node], type_ignores=[]),
                 str(EXPORT_MENU), "exec"), namespace)

    calls = []
    fake_menu = types.SimpleNamespace(launch_export=calls.append)
    monkeypatch.setitem(sys.modules, "cubit_export_menu", fake_menu)
    original_path = list(sys.path)
    try:
        for fmt in ("netgen", "gmsh", "nastran", "vtk", "femeem", "meg"):
            exec(namespace["_activate_code"](fmt), {})
    finally:
        sys.path[:] = original_path

    assert calls == ["netgen", "gmsh", "nastran", "vtk", "femeem", "meg"]
