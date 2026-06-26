"""Regression tests for the IH Run button after the workpiece .vol changes.

The original 2026-05-08 regression was reported on the standalone
``radia-ih`` panel: non-vacuum methods stayed disabled after the user
picked a valid workpiece ``.vol``.  The panel has since moved from the
legacy top-level ``AnalysisWindow._vol_edit`` row to a panel-owned
``IHPanel.wp_vol`` browse row under a working-folder root.  These tests
pin the current contract:

* the workpiece browse row re-inspects labels on text change;
* saved ``working_folder`` + relative ``panel.wp_vol`` restore together;
* constructor-supplied ``vol_path`` still overrides saved panel state.
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest


REPO = Path(__file__).resolve().parents[2]
SAMPLE_VOL = REPO / "tests" / "panels" / "test_3d_sibc_copper.vol"


@pytest.fixture
def isolated_settings(tmp_path, monkeypatch):
    """Redirect panel settings to a temporary directory."""
    import radia_gui_base
    monkeypatch.setattr(radia_gui_base, "_SETTINGS_DIR", str(tmp_path))
    return tmp_path


@pytest.fixture
def ih_window_clean(qapp, isolated_settings):
    """IHWindow with isolated settings and no constructor .vol."""
    from radia_ih import IHWindow
    win = IHWindow(vol_path="")
    yield win
    win.deleteLater()


def _set_wp_vol(win, path: Path) -> None:
    """Simulate the workpiece browse row accepting a path."""
    le = win._panel._widgets["wp_vol"]
    le.setText(win.display_path(str(path)))
    win._panel._on_wp_vol_changed_text(le.text())


class TestRunButtonAfterWorkpieceBrowse:
    """The panel-owned wp_vol row must drive label inspection."""

    @pytest.mark.parametrize("method_const", [
        "METHOD_PEEC_BEM",
        "METHOD_BEMA_BEM",
        "METHOD_PEEC_FEM_KELVIN",
    ])
    def test_run_disabled_on_fresh_launch_with_empty_vol(
            self, ih_window_clean, method_const):
        import radia_ih
        method = getattr(radia_ih, method_const)
        ih_window_clean._panel._method_combo.setCurrentText(method)
        assert ih_window_clean._panel._vol_mats is None
        assert ih_window_clean._panel.is_runnable() is False
        assert ih_window_clean._run_btn.isEnabled() is False

    @pytest.mark.parametrize("method_const", [
        "METHOD_PEEC_BEM",
        "METHOD_BEMA_BEM",
        "METHOD_PEEC_FEM_KELVIN",
        "METHOD_FEM_FULL",
    ])
    def test_browse_to_valid_vol_enables_run(
            self, ih_window_clean, method_const):
        import radia_ih
        method = getattr(radia_ih, method_const)
        ih_window_clean._panel._method_combo.setCurrentText(method)

        _set_wp_vol(ih_window_clean, SAMPLE_VOL)

        assert ih_window_clean._panel._vol_mats == {"coil", "kelvin", "air"}
        assert ih_window_clean._panel._vol_bnds == {
            "sibc", "source", "sink", "kelvin_int", "kelvin_ext", "default",
        }
        assert ih_window_clean._panel.is_runnable() is True
        assert ih_window_clean._run_btn.isEnabled() is True

    def test_browse_to_nonexistent_vol_keeps_run_disabled(
            self, ih_window_clean, tmp_path):
        import radia_ih
        ih_window_clean._panel._method_combo.setCurrentText(
            radia_ih.METHOD_PEEC_BEM)

        _set_wp_vol(ih_window_clean, SAMPLE_VOL)
        assert ih_window_clean._run_btn.isEnabled() is True

        bogus = tmp_path / "does_not_exist.vol"
        _set_wp_vol(ih_window_clean, bogus)
        assert ih_window_clean._panel._vol_mats is None
        assert ih_window_clean._panel._vol_bnds is None
        assert ih_window_clean._run_btn.isEnabled() is False


class TestSettingsRoundTrip:
    """Saved working folder and relative browse rows restore together."""

    def test_init_uses_restored_relative_wp_vol(
            self, qapp, isolated_settings):
        from radia_ih import IHWindow, METHOD_PEEC_BEM

        json_path = Path(isolated_settings) / "radia_ih.json"
        json_path.write_text(json.dumps({
            "working_folder": str(SAMPLE_VOL.parent),
            "vol": "",
            "panel": {
                "method": METHOD_PEEC_BEM,
                "wp_vol": SAMPLE_VOL.name,
            },
        }), encoding="utf-8")

        win = IHWindow(vol_path="")
        try:
            assert Path(win.working_folder) == SAMPLE_VOL.parent
            assert win._panel._widgets["wp_vol"].text() == SAMPLE_VOL.name
            assert Path(win._panel.wp_vol_path()) == SAMPLE_VOL
            assert win._panel._vol_mats == {"coil", "kelvin", "air"}
            assert win._run_btn.isEnabled() is True
        finally:
            win.deleteLater()

    def test_save_browse_restore_round_trip(
            self, qapp, isolated_settings):
        from radia_ih import IHWindow, METHOD_PEEC_BEM

        win1 = IHWindow(vol_path="")
        try:
            win1._folder_edit.setText(str(SAMPLE_VOL.parent))
            win1._on_folder_changed(str(SAMPLE_VOL.parent))
            win1._panel._method_combo.setCurrentText(METHOD_PEEC_BEM)
            _set_wp_vol(win1, SAMPLE_VOL)
            assert win1._panel._widgets["wp_vol"].text() == SAMPLE_VOL.name
            assert win1._run_btn.isEnabled() is True
            win1._save_settings()
        finally:
            win1.deleteLater()

        win2 = IHWindow(vol_path="")
        try:
            assert Path(win2.working_folder) == SAMPLE_VOL.parent
            assert win2._panel._widgets["wp_vol"].text() == SAMPLE_VOL.name
            assert win2._panel._method_combo.currentText() == METHOD_PEEC_BEM
            assert win2._run_btn.isEnabled() is True
        finally:
            win2.deleteLater()

    def test_constructor_vol_overrides_restored_wp_vol(
            self, qapp, isolated_settings, tmp_path):
        from radia_ih import IHWindow, METHOD_PEEC_BEM

        stale = tmp_path / "stale.vol"
        json_path = Path(isolated_settings) / "radia_ih.json"
        json_path.write_text(json.dumps({
            "working_folder": str(tmp_path),
            "vol": "",
            "panel": {
                "method": METHOD_PEEC_BEM,
                "wp_vol": str(stale),
            },
        }), encoding="utf-8")

        win = IHWindow(vol_path=str(SAMPLE_VOL))
        try:
            assert Path(win.working_folder) == SAMPLE_VOL.parent
            assert Path(win._panel.wp_vol_path()) == SAMPLE_VOL
            assert win._run_btn.isEnabled() is True
        finally:
            win.deleteLater()


class TestWorkpieceBrowseContract:
    """The current browse hook lives on IHPanel, not AnalysisWindow."""

    def test_analysis_window_folder_change_refreshes_browse_rows(
            self, ih_window_clean):
        le = ih_window_clean._panel._widgets["wp_vol"]
        le.setText(str(SAMPLE_VOL))

        ih_window_clean._folder_edit.setText(str(SAMPLE_VOL.parent))
        ih_window_clean._on_folder_changed(str(SAMPLE_VOL.parent))

        assert le.text() == SAMPLE_VOL.name
        assert Path(ih_window_clean._panel.wp_vol_path()) == SAMPLE_VOL

    def test_workpiece_browse_row_reloads_labels(
            self, ih_window_clean, monkeypatch):
        from PySide6.QtWidgets import QFileDialog
        import radia_ih

        monkeypatch.setattr(QFileDialog, "getOpenFileName",
                            staticmethod(lambda *a, **kw:
                                         (str(SAMPLE_VOL), "")))
        ih_window_clean._panel._method_combo.setCurrentText(
            radia_ih.METHOD_PEEC_BEM)
        le = ih_window_clean._panel._widgets["wp_vol"]

        ih_window_clean._panel._do_browse(
            le, le.property("_radia_browse_filter"))

        assert Path(ih_window_clean._panel.wp_vol_path()) == SAMPLE_VOL
        assert ih_window_clean._panel._vol_mats == {"coil", "kelvin", "air"}
        assert ih_window_clean._run_btn.isEnabled() is True
