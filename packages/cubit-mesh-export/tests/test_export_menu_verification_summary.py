"""The toolbar's Netgen export must leave the check-vol verdict in the console log."""

import importlib
import sys
import types


def _load_menu(monkeypatch):
    # The menu runs inside Cubit's PySide6; the summary helper needs no Qt.
    widgets = types.ModuleType("PySide6.QtWidgets")
    for name in ("QApplication", "QCheckBox", "QComboBox", "QDialog",
                 "QDialogButtonBox", "QFileDialog", "QFormLayout", "QGroupBox",
                 "QHBoxLayout", "QLabel", "QLineEdit", "QMainWindow",
                 "QMessageBox", "QPushButton", "QTableWidget",
                 "QTableWidgetItem", "QVBoxLayout"):
        setattr(widgets, name, type(name, (), {}))
    core = types.ModuleType("PySide6.QtCore")
    core.Qt = type("Qt", (), {})
    package = types.ModuleType("PySide6")
    monkeypatch.setitem(sys.modules, "PySide6", package)
    monkeypatch.setitem(sys.modules, "PySide6.QtWidgets", widgets)
    monkeypatch.setitem(sys.modules, "PySide6.QtCore", core)
    monkeypatch.delitem(
        sys.modules, "cubit_mesh_export.cubit_gui.cubit_export_menu", raising=False)
    return importlib.import_module("cubit_mesh_export.cubit_gui.cubit_export_menu")


def test_failed_report_prints_failure_and_its_warnings(monkeypatch):
    menu = _load_menu(monkeypatch)
    text = menu._verification_summary(
        1, {"passed": False, "n_elements": 0, "n_points": 1331, "order": 3,
            "warnings": ["no volume elements"]}, "o:/m.vol.check.json")
    assert text.splitlines()[0] == (
        "Verification FAILED: 0 elements, 1331 points, order 3 "
        "(report: o:/m.vol.check.json)")
    assert "  warning: no volume elements" in text.splitlines()


def test_passed_report_needs_zero_exit_code(monkeypatch):
    menu = _load_menu(monkeypatch)
    report = {"passed": True, "n_elements": 1000, "n_points": 6131, "order": 3}
    assert menu._verification_summary(0, report, "r").startswith(
        "Verification passed: 1000 elements")
    assert menu._verification_summary(2, report, "r").startswith(
        "Verification FAILED")


def test_error_report_prints_the_error(monkeypatch):
    menu = _load_menu(monkeypatch)
    text = menu._verification_summary(
        2, {"passed": False, "error": "cannot read mesh"}, "r")
    assert "  error: cannot read mesh" in text.splitlines()
    assert "? elements" in text
