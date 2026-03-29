"""
GUI test: InductanceDialog operated programmatically (no Cubit, no human).

Uses PySide6 QTest to simulate button clicks, combo selection, text input.
cubit module is mocked so the dialog runs standalone.

Usage:
    python tests/panels/test_panel_gui.py          # run all
    python tests/panels/test_panel_gui.py -v        # verbose
    pytest tests/panels/test_panel_gui.py -v        # via pytest
"""

import os
import sys
import types
import json
import pytest

# ============================================================
# Mock cubit module (must be done before importing register_toolbar)
# ============================================================
def _install_cubit_mock(block_names=None):
    """Install a mock cubit module into sys.modules.

    Args:
        block_names: dict of {name: block_id} to simulate Cubit blocks.
                     e.g. {"source": 1, "sink": 2, "workpiece": 3}
    """
    if block_names is None:
        block_names = {}

    mock = types.ModuleType("cubit")

    # Reverse map: id -> name
    _id_to_name = {v: k for k, v in block_names.items()}

    mock.get_block_id_list = lambda: tuple(block_names.values())
    mock.get_exodus_entity_name = lambda kind, bid: _id_to_name.get(bid, "")
    mock.get_entities = lambda kind: ()
    mock.cmd = lambda s: None
    mock.get_file_name = lambda: ""
    mock.get_block_volumes = lambda bid: ()

    # Surface mock for workpiece
    class _MockSurface:
        def center_point(self): return (0, 0, 0)
        def normal_at(self, pt): return (0, 0, 1)
        def area(self): return 1e-4
        def is_planar(self): return False

    mock.surface = lambda sid: _MockSurface()
    mock.get_relatives = lambda *a: ()
    mock.get_volume_tets = lambda vid: (1, 2, 3)
    mock.get_volume_hexes = lambda vid: ()

    sys.modules["cubit"] = mock
    return mock


# ============================================================
# Fixture: QApplication (one per session)
# ============================================================
@pytest.fixture(scope="session")
def qapp():
    from PySide6.QtWidgets import QApplication
    app = QApplication.instance()
    if app is None:
        app = QApplication(sys.argv)
    return app


@pytest.fixture
def dialog_no_workpiece(qapp):
    """InductanceDialog with source+sink but no workpiece."""
    _install_cubit_mock({"source": 1, "sink": 2, "conductor": 3})
    # Reload to pick up mock
    import importlib
    _repo = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    panels_dir = os.path.join(_repo, "src", "radia", "panels")
    if panels_dir not in sys.path:
        sys.path.insert(0, panels_dir)
    pkg_root = os.path.join(_repo, "src", "radia")
    if pkg_root not in sys.path:
        sys.path.insert(0, pkg_root)

    if "register_toolbar" in sys.modules:
        del sys.modules["register_toolbar"]
    import register_toolbar
    dlg = register_toolbar.InductanceDialog()
    yield dlg
    dlg.close()


@pytest.fixture
def dialog_with_workpiece(qapp):
    """InductanceDialog with source+sink+workpiece."""
    _install_cubit_mock({"source": 1, "sink": 2, "conductor": 3, "workpiece": 4})
    import importlib
    _repo = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    panels_dir = os.path.join(_repo, "src", "radia", "panels")
    if panels_dir not in sys.path:
        sys.path.insert(0, panels_dir)
    pkg_root = os.path.join(_repo, "src", "radia")
    if pkg_root not in sys.path:
        sys.path.insert(0, pkg_root)

    if "register_toolbar" in sys.modules:
        del sys.modules["register_toolbar"]
    import register_toolbar
    dlg = register_toolbar.InductanceDialog()
    yield dlg
    dlg.close()


# ============================================================
# Tests: Dialog construction and block detection
# ============================================================
class TestDialogConstruction:

    def test_dialog_opens(self, dialog_no_workpiece):
        """Dialog should instantiate without error."""
        dlg = dialog_no_workpiece
        assert dlg.windowTitle() is not None
        assert dlg.isVisible() is False  # not shown yet

    def test_source_sink_detected(self, dialog_no_workpiece):
        """Source and sink labels should be green."""
        dlg = dialog_no_workpiece
        assert dlg.source_label.text() == "source"
        assert dlg.sink_label.text() == "sink"
        assert "green" in dlg.source_label.styleSheet()
        assert "green" in dlg.sink_label.styleSheet()

    def test_workpiece_not_found(self, dialog_no_workpiece):
        """Without workpiece block, ESIM group should be hidden."""
        dlg = dialog_no_workpiece
        assert dlg.workpiece_label.text() == "(none)"
        # isVisible() depends on parent visibility; use isHidden() for explicit hide
        assert dlg.esim_group.isHidden() is True

    def test_workpiece_found(self, dialog_with_workpiece):
        """With workpiece block, ESIM group should not be hidden."""
        dlg = dialog_with_workpiece
        assert dlg.workpiece_label.text() == "workpiece"
        assert "green" in dlg.workpiece_label.styleSheet()
        # setVisible(True) was called; isHidden() checks the explicit flag
        assert dlg.esim_group.isHidden() is False


# ============================================================
# Tests: Widget interaction (QTest-style programmatic clicks)
# ============================================================
class TestWidgetInteraction:

    def test_sigma_editable(self, dialog_with_workpiece):
        """Sigma field should accept text input."""
        dlg = dialog_with_workpiece
        dlg.sigma_edit.setText("5.8e7")
        assert dlg.sigma_edit.text() == "5.8e7"

    def test_model_combo_options(self, dialog_with_workpiece):
        """Model combo should have ESIM and Dowell."""
        dlg = dialog_with_workpiece
        items = [dlg.model_combo.itemText(i)
                 for i in range(dlg.model_combo.count())]
        assert "ESIM" in items
        assert "Dowell" in items

    def test_frequency_editable(self, dialog_with_workpiece):
        """Frequency field should accept text input."""
        dlg = dialog_with_workpiece
        dlg.freq_edit.setText("100000")
        assert dlg.freq_edit.text() == "100000"

    def test_curve_spin_range(self, dialog_with_workpiece):
        """Curve order spin should be 1-2."""
        dlg = dialog_with_workpiece
        assert dlg.curve_spin.minimum() == 1
        assert dlg.curve_spin.maximum() == 2

    def test_fes_spin_range(self, dialog_with_workpiece):
        """FES order spin should be 0-2."""
        dlg = dialog_with_workpiece
        assert dlg.fes_spin.minimum() == 0
        assert dlg.fes_spin.maximum() == 2


# ============================================================
# Tests: Result display
# ============================================================
class TestResultDisplay:

    def test_display_inductance_only(self, dialog_no_workpiece):
        """Display result without ESIM data."""
        dlg = dialog_no_workpiece
        data = {
            "inductance_H": 150e-9,
            "n_dofs": 5000,
            "n_faces": 3000,
            "t_export": 0.5,
            "t_assembly": 10.0,
            "t_lu": 5.0,
            "source_area": 7.8e-5,
            "sink_area": 7.8e-5,
            "surface_area": 3.1e-3,
            "constraint_residual": 1e-14,
            "curve_order": 2,
            "fes_order": 0,
        }
        dlg._display_result(data)

        # Check table has rows
        assert dlg.result_table.rowCount() >= 10
        # First row should be inductance
        assert "Inductance" in dlg.result_table.item(0, 0).text()
        assert "nH" in dlg.result_table.item(0, 1).text()

    def test_display_with_esim(self, dialog_with_workpiece):
        """Display result with ESIM workpiece data."""
        dlg = dialog_with_workpiece
        data = {
            "inductance_H": 150e-9,
            "n_dofs": 5000,
            "n_faces": 3000,
            "t_export": 0.5,
            "t_assembly": 10.0,
            "t_lu": 5.0,
            "source_area": 7.8e-5,
            "sink_area": 7.8e-5,
            "surface_area": 3.1e-3,
            "constraint_residual": 1e-14,
            "curve_order": 2,
            "fes_order": 0,
            "wp_R_effective": 6.2e-4,
            "wp_X_effective": 1.4e-3,
            "wp_P_total": 3.1e-4,
            "wp_Q_total": 7.2e-4,
            "wp_frequency": 50000,
            "wp_delta_min": 0.04e-3,
            "wp_delta_max": 0.04e-3,
            "wp_n_panels": 128,
            "wp_model": "esim",
            "wp_material": "steel",
        }
        dlg._display_result(data)

        # Should have more rows than inductance-only
        n_rows = dlg.result_table.rowCount()
        assert n_rows >= 20, f"Expected >= 20 rows with ESIM, got {n_rows}"

        # Find R (workpiece) row
        found_R = False
        found_total_Z = False
        for i in range(n_rows):
            item = dlg.result_table.item(i, 0)
            if item and "R (workpiece)" in item.text():
                found_R = True
                val = dlg.result_table.item(i, 1).text()
                assert "Ohm" in val
            if item and "|Z_total|" in item.text():
                found_total_Z = True

        assert found_R, "R (workpiece) row not found"
        assert found_total_Z, "|Z_total| row not found"


# ============================================================
# Tests: Solve button state
# ============================================================
class TestSolveButtonState:

    def test_solve_L_enabled_with_blocks(self, dialog_no_workpiece):
        """Solve L button should reflect source+sink+python availability."""
        dlg = dialog_no_workpiece
        if dlg._ext_python:
            assert dlg.solve_btn.isEnabled() is True
        else:
            assert dlg.solve_btn.isEnabled() is False

    def test_esim_group_hidden_without_workpiece(self, dialog_no_workpiece):
        """ESIM group hidden without workpiece block."""
        dlg = dialog_no_workpiece
        assert dlg.esim_group.isHidden() is True

    def test_esim_group_visible_with_workpiece(self, dialog_with_workpiece):
        """ESIM group visible when workpiece block found."""
        dlg = dialog_with_workpiece
        assert dlg.esim_group.isHidden() is False

    def test_solve_disabled_without_blocks(self, qapp):
        """No blocks -> solve disabled."""
        _install_cubit_mock({})  # no blocks
        _repo = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
        panels_dir = os.path.join(_repo, "src", "radia", "panels")
        if panels_dir not in sys.path:
            sys.path.insert(0, panels_dir)
        pkg_root = os.path.join(_repo, "src", "radia")
        if pkg_root not in sys.path:
            sys.path.insert(0, pkg_root)
        if "register_toolbar" in sys.modules:
            del sys.modules["register_toolbar"]
        import register_toolbar
        dlg = register_toolbar.InductanceDialog()
        assert dlg.solve_btn.isEnabled() is False
        dlg.close()


# ============================================================
if __name__ == "__main__":
    pytest.main([__file__, "-v", "--tb=short"])
