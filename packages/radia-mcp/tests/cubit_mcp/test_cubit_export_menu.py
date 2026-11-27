"""Guards: the Netgen Export menu verifies with radia.vol_check, routes Kelvin through
radia_mcp.cubit.kelvin, leaves the verdict in the console log and loads no radia_mcp code.

Assumes: no Cubit and no Qt; PySide6 is replaced by empty stand-in classes, so
only the Qt-free helpers of the menu module are exercised.
"""
from __future__ import annotations

import importlib.util
import subprocess
import sys
import types
from pathlib import Path

import pytest

import radia_mcp.cubit
from radia_mcp.cubit import kelvin

GUI = Path(radia_mcp.cubit.__file__).resolve().parent / "gui"
MENU = GUI / "cubit_export_menu.py"


def _fake_pyside(monkeypatch):
    widgets = types.ModuleType("PySide6.QtWidgets")
    for name in ("QApplication", "QCheckBox", "QComboBox", "QDialog",
                 "QDialogButtonBox", "QFileDialog", "QFormLayout", "QGroupBox",
                 "QHBoxLayout", "QLabel", "QLineEdit", "QMainWindow",
                 "QMessageBox", "QPushButton", "QTableWidget",
                 "QTableWidgetItem", "QVBoxLayout"):
        setattr(widgets, name, type(name, (), {}))
    core = types.ModuleType("PySide6.QtCore")
    core.Qt = type("Qt", (), {})
    core.QTimer = type("QTimer", (), {})
    gui = types.ModuleType("PySide6.QtGui")
    gui.QAction = type("QAction", (), {})
    monkeypatch.setitem(sys.modules, "PySide6", types.ModuleType("PySide6"))
    monkeypatch.setitem(sys.modules, "PySide6.QtWidgets", widgets)
    monkeypatch.setitem(sys.modules, "PySide6.QtCore", core)
    monkeypatch.setitem(sys.modules, "PySide6.QtGui", gui)


@pytest.fixture
def menu(monkeypatch):
    """The menu loaded by path, as register_toolbar.py and the toolbar script do."""
    _fake_pyside(monkeypatch)
    spec = importlib.util.spec_from_file_location("cubit_export_menu", MENU)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_failed_report_prints_failure_and_its_warnings(menu):
    text = menu._verification_summary(
        1, {"passed": False, "n_elements": 0, "n_points": 1331, "order": 3,
            "warnings": ["no volume elements"]}, "o:/m.vol.check.json")
    assert text.splitlines()[0] == (
        "Verification FAILED: 0 elements, 1331 points, order 3 "
        "(report: o:/m.vol.check.json)")
    assert "  warning: no volume elements" in text.splitlines()


def test_passed_report_needs_zero_exit_code(menu):
    report = {"passed": True, "n_elements": 1000, "n_points": 6131, "order": 3}
    assert menu._verification_summary(0, report, "r").startswith(
        "Verification passed: 1000 elements")
    assert menu._verification_summary(2, report, "r").startswith(
        "Verification FAILED")


def test_error_report_prints_the_error(menu):
    text = menu._verification_summary(
        2, {"passed": False, "error": "cannot read mesh"}, "r")
    assert "  error: cannot read mesh" in text.splitlines()
    assert "? elements" in text


def test_verification_runs_radia_vol_check_isolated(menu):
    command = menu.verification_command("C:/py/python.exe", "o:/m.vol", "o:/m.vol.check.json")
    assert command == ["C:/py/python.exe", "-I", "-B", "-m", "radia.vol_check", "o:/m.vol",
                       "--format", "json", "--report-json", "o:/m.vol.check.json"], command


def test_export_command_carries_no_kelvin_keywords(menu):
    # The native plugin has no add_kelvin keyword: Kelvin runs before export.
    assert menu.export_command("o:/m.vol", 3) == 'export netgen "o:/m.vol" order 3 overwrite'


@pytest.mark.parametrize("environment", [None, "python", "Z:/no/such/python.exe"])
def test_external_python_must_be_an_existing_absolute_path(menu, monkeypatch, environment):
    if environment is None:
        monkeypatch.delenv(menu.PYTHON_ENV, raising=False)
    else:
        monkeypatch.setenv(menu.PYTHON_ENV, environment)
    with pytest.raises(RuntimeError, match=menu.PYTHON_ENV):
        menu._external_python()


def test_external_python_accepts_the_recorded_interpreter(menu, monkeypatch):
    assert menu.PYTHON_ENV == "RADIA_MCP_CUBIT_PYTHON"
    monkeypatch.setenv(menu.PYTHON_ENV, sys.executable)
    assert menu._external_python() == sys.executable


def test_dialog_kelvin_fields_map_onto_the_kelvin_contract(menu):
    options = menu.kelvin_options(True, " 0.03 ", "bn", "off", "ht")
    normalized = kelvin.validate_options(options)
    assert normalized["kelvin_mesh"] == 0.03, normalized
    assert normalized["reduction"] == {"x": "bn=0", "z": "ht=0"}, normalized
    assert kelvin.validate_options(menu.kelvin_options(True, "", "off", "off", "off"))[
        "kelvin_mesh"] is None
    with pytest.raises(ValueError):
        kelvin.validate_options(menu.kelvin_options(True, "abc", "off", "off", "off"))


def test_menu_loads_the_package_kelvin_module_by_path(menu):
    loaded = menu._kelvin_module()
    assert Path(loaded.__file__).resolve() == Path(kelvin.__file__).resolve()
    assert loaded.__name__ != kelvin.__name__  # a by-path copy, not radia_mcp.cubit.kelvin


def test_gui_import_loads_no_radia_mcp_code_beyond_its_packages():
    probe = (
        "import sys, types\n"
        "for n in ('PySide6','PySide6.QtWidgets','PySide6.QtCore','PySide6.QtGui'):\n"
        "    sys.modules[n] = types.ModuleType(n)\n"
        "w = sys.modules['PySide6.QtWidgets']\n"
        "for n in ('QApplication','QCheckBox','QComboBox','QDialog','QDialogButtonBox',"
        "'QFileDialog','QFormLayout','QGroupBox','QHBoxLayout','QLabel','QLineEdit',"
        "'QMainWindow','QMessageBox','QPushButton','QTableWidget','QTableWidgetItem','QVBoxLayout'):\n"
        "    setattr(w, n, type(n, (), {}))\n"
        "sys.modules['PySide6.QtCore'].Qt = type('Qt', (), {})\n"
        "sys.modules['PySide6.QtCore'].QTimer = type('QTimer', (), {})\n"
        "sys.modules['PySide6.QtGui'].QAction = type('QAction', (), {})\n"
        "import radia_mcp.cubit.gui.cubit_export_menu\n"
        "print(sorted(m for m in sys.modules if m.startswith(('radia', 'ngsolve', 'netgen', 'numpy'))))\n"
    )
    run = subprocess.run([sys.executable, "-c", probe], capture_output=True, text=True, timeout=120)
    assert run.returncode == 0, run.stderr
    loaded = eval(run.stdout.strip().splitlines()[-1])
    assert loaded == ["radia_mcp", "radia_mcp.cubit", "radia_mcp.cubit.gui",
                      "radia_mcp.cubit.gui.cubit_export_menu"], loaded


def test_menu_owns_one_action_and_never_walks_cubit_menus():
    source = MENU.read_text(encoding="utf-8")
    assert "app.aboutToQuit.connect(remove_menu)" in source
    assert "bar.addAction(action)" in source
    for forbidden in ("findChildren", ".actions()", "createPopupMenu", "import radia_mcp",
                      "from radia_mcp"):
        assert forbidden not in source, forbidden
